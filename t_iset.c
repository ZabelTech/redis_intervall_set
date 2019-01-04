#include "redismodule.h"
#include "dict.h"

#include <math.h>
#include <string.h>

//=============================

static RedisModuleType *ISet;

int dictEncObjKeyCompare(void *privdata, const void *key1, const void *key2) {
  DICT_NOTUSED(privdata);
  return RedisModule_StringCompare((RedisModuleString *) key1,
				   (RedisModuleString *) key2) == 0;
}

uint64_t dictEncObjHash(const void *key) {
  size_t len;
  const char * keyStr;
  keyStr = RedisModule_StringPtrLen((RedisModuleString *)key,&len);
  return dictGenHashFunction(keyStr, len);
}

void dictObjectDestructor(void *privdata, void *val) {
    DICT_NOTUSED(privdata);
    if (val == NULL) return;
    RedisModule_Free((RedisModuleString *) val);
}

dictType isetDictType = {
    dictEncObjHash,            /* hash function */
    NULL,                      /* key dup */
    NULL,                      /* val dup */
    dictEncObjKeyCompare,      /* key compare */
    dictObjectDestructor,      /* key destructor */
    NULL                       /* val destructor */
};
//=============================

/*-----------------------------------------------------------------------------
 * Interval set API
 *----------------------------------------------------------------------------*/

/* ISETs are sets using two data structures to hold the same elements
 * in order to get O(log(N)) INSERT and REMOVE operations into an interval
 * range data structure.
 *
 * The elements are added to an hash table mapping Redis objects to intervals.
 * At the same time the elements are added to an augmented AVL tree that maps
 * intervals to Redis objects. */

typedef struct avlNode {
       RedisModuleString *obj;
       double scores[2];
       double subLeftMax, subRightMax;
       char balance;
       struct avlNode *left, *right, *parent, *next;
} avlNode;

typedef struct avl {
       struct avlNode *root;
       dict *dict;
       unsigned long size;
} avl;

avl *avlCreate(void) {
    avl *avltree;

    avltree = RedisModule_Alloc(sizeof(*avltree));
    avltree->size = 0;
    avltree->root = NULL;

    avltree->dict = dictCreate(&isetDictType,NULL);

    return avltree;
}

avlNode *avlCreateNode(RedisModuleCtx *ctx, double lscore, double rscore, RedisModuleString *obj) {
    avlNode *an = RedisModule_Alloc(sizeof(*an));
    an->scores[0] = lscore;
    an->scores[1] = rscore;
    an->subLeftMax = -INFINITY;
    an->subRightMax = -INFINITY;
    an->balance = 0;
    an->left = NULL;
    an->right = NULL;
    an->parent = NULL;
    an->next = NULL;
    an->obj = obj;

    if (obj)
	RedisModule_RetainString(ctx, obj);

    return an;
}

void avlFreeNode(RedisModuleCtx *ctx, avlNode *node, int eraseList) {
    if (eraseList && node->next)
	avlFreeNode(ctx, node->next, eraseList);
    if (node->obj)
	RedisModule_FreeString(ctx,node->obj);
    if (node->left)
	avlFreeNode(ctx, node->left, eraseList);
    if (node->right)
	avlFreeNode(ctx, node->right, eraseList);
    RedisModule_Free(node);
}

void avlFreeTree(avlNode *node, int eraseList) {
    if (eraseList && node->next)
	avlFreeTree(node->next, eraseList);
    if (node->obj)
	RedisModule_Free(node->obj);
    if (node->left)
	avlFreeTree(node->left, eraseList);
    if (node->right)
	avlFreeTree(node->right, eraseList);
    RedisModule_Free(node);
}

void avlFree(void *o) {
    avl *tree = o;
    if (tree->root != NULL)
	avlFreeTree(tree->root,1);
    dictRelease(tree->dict);
    RedisModule_Free(tree);
}

int avlNodeCmp(avlNode *a, avlNode *b) {
    if (a->scores[0] < b->scores[0])
	return -1;
    else if (a->scores[0] > b->scores[0])
	return 1;
    else {
	if (a->scores[1] > b->scores[1])
	    return -1;
	else if (a->scores[1] < b->scores[1])
	    return 1;
	else
	    return 0;
    }
}

void avlLeftRotation(avl * tree, avlNode *locNode) {
    avlNode *newRoot = locNode->right;
    locNode->right = newRoot->left;
    if (locNode->right) locNode->right->parent = locNode;
    newRoot->left = locNode;

    newRoot->parent = locNode->parent;
    locNode->parent = newRoot;
    if (newRoot->parent) {
	if (avlNodeCmp(newRoot->parent,newRoot) > -1)
	    newRoot->parent->left = newRoot;
	else
	    newRoot->parent->right = newRoot;
    }
    // New root
    else {
	tree->root = newRoot;
    }
}

void avlRightRotation(avl * tree, avlNode *locNode) {
    avlNode *newRoot = locNode->left;
    locNode->left = newRoot->right;
    if (locNode->left) locNode->left->parent = locNode;
    newRoot->right = locNode;

    newRoot->parent = locNode->parent;
    locNode->parent = newRoot;
    if (newRoot->parent) {
	if(avlNodeCmp(newRoot->parent,newRoot) > -1)
	    newRoot->parent->left = newRoot;
	else
	    newRoot->parent->right = newRoot;
    }
    // New root
    else {
	tree->root = newRoot;
    }
}

void avlResetBalance(avlNode *locNode) {
    switch(locNode->balance) {
	case -1:
	locNode->left->balance = 0;
	locNode->right->balance = 1;
	break;
	case 0:
	locNode->left->balance = 0;
	locNode->right->balance = 0;
	break;
	case 1:
	locNode->left->balance = -1;
	locNode->right->balance = 0;
	break;
    }
    locNode->balance = 0;
}

void avlUpdateMaxScores(avlNode *locNode) {
    double oldNodeMax;

    while (locNode) {
	if (locNode->left) {
	    oldNodeMax = locNode->left->scores[1];
	    oldNodeMax = (oldNodeMax > locNode->left->subLeftMax) ? oldNodeMax : locNode->left->subLeftMax;
	    oldNodeMax = (oldNodeMax > locNode->left->subRightMax) ? oldNodeMax : locNode->left->subRightMax;
	    locNode->subLeftMax = oldNodeMax;
	}
	else {
	    locNode->subLeftMax = -INFINITY;
	}
	if (locNode->right) {
	    oldNodeMax = locNode->right->scores[1];
	    oldNodeMax = (oldNodeMax > locNode->right->subLeftMax) ? oldNodeMax : locNode->right->subLeftMax;
	    oldNodeMax = (oldNodeMax > locNode->right->subRightMax) ? oldNodeMax : locNode->right->subRightMax;
	    locNode->subRightMax = oldNodeMax;
	}
	else {
	    locNode->subRightMax = -INFINITY;
	}
	locNode = locNode->parent;
    }
}

int avlInsertNode(avl * tree, avlNode *locNode, avlNode *insertNode) {
    int diff = avlNodeCmp(locNode, insertNode);
    /* Insert in the left node */
    if (diff > 0) {
	if (!locNode->left) {
	    locNode->left = insertNode;
	    insertNode->parent = locNode;
	    locNode->balance = locNode->balance - 1;
	    avlUpdateMaxScores(locNode);
	    return locNode->balance ? 1 : 0;
	}
	else {
	    // Left node is occupied, insert it into the subtree
	    if (avlInsertNode(tree,locNode->left,insertNode)) {
		locNode->balance = locNode->balance - 1;
		if (locNode->balance == 0)
		    return 0;
		else if (locNode->balance == -1)
		    return 1;

		// Tree is unbalanced at this point
		// Case 1 at http://www.stanford.edu/~blp/avl/libavl.html/Rebalancing-AVL-Trees.html#index-rebalance-after-AVL-insertion-236
		if (locNode->left->balance < 0) {
		    // Left-Left, single right rotation needed
		    avlRightRotation(tree,locNode);

		    //Both locNode and its parent have a zero balance; see note at the link above
		    locNode->balance = 0;
		    locNode->parent->balance = 0;

		    locNode->subLeftMax = locNode->parent->subRightMax;
		    //XXX: What is this, I don't even?
		    //     locNode->parent->subRightMax should be the max of the subtree
		    //     *rooted at locNode*, right?
		    locNode->parent->subRightMax = -INFINITY;
		}
		else {
		    // Left-Right, left rotation then right rotation needed
		    avlLeftRotation(tree,locNode->left);
		    avlRightRotation(tree,locNode);
		    avlResetBalance(locNode->parent);

		    locNode->subLeftMax = locNode->parent->subRightMax;
		    locNode->parent->left->subRightMax = locNode->parent->subLeftMax;
		    locNode->parent->subRightMax = -INFINITY;
		    locNode->parent->subLeftMax = -INFINITY;
		}

		avlUpdateMaxScores(locNode->parent);
	    }
	    return 0;
	}
    }
    /* Insert in the right node */
    else if (diff < 0) {
	if (!locNode->right) {
	    locNode->right = insertNode;
	    insertNode->parent = locNode;
	    locNode->balance = locNode->balance + 1;
	    avlUpdateMaxScores(locNode);
	    return locNode->balance ? 1 : 0;
	}
	else {
	    // Right node is occupied, insert it into the subtree
	    if (avlInsertNode(tree,locNode->right,insertNode)) {
		locNode->balance = locNode->balance + 1;
		if (locNode->balance == 0)
		    return 0;
		else if (locNode->balance == 1)
		    return 1;

		// Tree is unbalanced at this point
		if (locNode->right->balance > 0) {
		    // Right-Right, single left rotation needed
		    avlLeftRotation(tree,locNode);

		    locNode->balance = 0;
		    locNode->parent->balance = 0;

		    locNode->subRightMax = locNode->parent->subLeftMax;
		    locNode->parent->subLeftMax = -INFINITY;
		}
		else {
		    // Right-Left, right rotation then left rotation needed
		    avlRightRotation(tree,locNode->right);
		    avlLeftRotation(tree,locNode);
		    avlResetBalance(locNode->parent);

		    locNode->subRightMax = locNode->parent->subLeftMax;
		    locNode->parent->right->subLeftMax = locNode->parent->subRightMax;
		    locNode->parent->subRightMax = -INFINITY;
		    locNode->parent->subLeftMax = -INFINITY;
		}

		avlUpdateMaxScores(locNode->parent);
	    }
	    return 0;
	}
    }
    // These nodes have the same range. We're going to assume that the RedisModuleString hasn't been
    // added before to this range, as the caller to avlInsert should check this
    else {
	avlNode * tail = locNode;
	while (tail->next != NULL)
	    tail = tail->next;
	tail->next = insertNode;
	return 0;
    }
}

avlNode *avlInsert(RedisModuleCtx *ctx, avl *tree, double lscore, double rscore, RedisModuleString *obj) {
    avlNode *an = avlCreateNode(ctx, lscore, rscore, obj);

    if (!tree->root)
	tree->root = an;
    else
	avlInsertNode(tree, tree->root, an);

    tree->size = tree->size + 1;

    return an;
}

void avlRemoveFromParent(avl * tree, avlNode *locNode, avlNode *replacementNode) {
    if (locNode->parent) {
	if (locNode->parent->left == locNode)
	    locNode->parent->left = replacementNode;
	else
	    locNode->parent->right = replacementNode;
    }
    else {
	tree->root = replacementNode;
    }
}

int avlRemoveNode(RedisModuleCtx *ctx, avl * tree, avlNode *locNode, avlNode *delNode, char freeNodeMem, int * removed) {
    int diff = avlNodeCmp(locNode, delNode);
    int heightDelta;
    avlNode *replacementNode;

    // This is the node we want removed
    if (diff == 0) {
	// First check to see if there are more than one element being stored here.
	// If so, find the element, remove it, and update the pointers appropriately.
	// If not, we can assume that this element is the one desired to be removed,
	// as the caller to avlRemoveNode should check the dict first to ensure the
	// obj exists at this point
	if (locNode->next && freeNodeMem) {
	    avlNode *removeNode = locNode;
	    avlNode *prevNode = NULL;

	    // Find the node where the node obj data matches the delNode obj data
	    while (RedisModule_StringCompare(removeNode->obj,delNode->obj) != 0) {
		prevNode = removeNode;
		removeNode = removeNode->next;
	    }
	    // If the node to be removed is the head, we need to update the locNode
	    if (removeNode == locNode) {
		locNode->next->parent = locNode->parent;
		locNode->next->left = locNode->left;
		locNode->next->right = locNode->right;
		locNode->next->balance = locNode->balance;
		locNode->next->subLeftMax = locNode->subLeftMax;
		locNode->next->subRightMax = locNode->subRightMax;

		// Update the parent and children
		avlRemoveFromParent(tree,locNode,locNode->next);
		if (locNode->left)
		    locNode->left->parent = locNode->next;
		if (locNode->right)
		    locNode->right->parent = locNode->next;

		locNode->right = NULL;
		locNode->left = NULL;
		avlFreeNode(ctx, locNode,0);
		*removed = 1;
		return 0;
	    }
	    else {
		prevNode->next = removeNode->next;
		avlFreeNode(ctx, removeNode,0);
		*removed = 1;
		return 0;
	    }
	}
	else {
	    // Remove if leaf node or replace with child if only one child
	    if (!locNode->left) {
		if (!locNode->right) {
		    avlRemoveFromParent(tree,locNode,NULL);
		    if (locNode->parent)
			avlUpdateMaxScores(locNode->parent);
		    if (freeNodeMem)
			avlFreeNode(ctx, locNode,0);
		    *removed = 1;
		    return -1;
		}
		avlRemoveFromParent(tree,locNode,locNode->right);
		locNode->right->parent = locNode->parent;
		if (locNode->parent)
		    avlUpdateMaxScores(locNode->parent);
		locNode->right = NULL;
		if (freeNodeMem)
		   avlFreeNode(ctx, locNode,0);
		*removed = 1;
		return -1;
	    }
	    if (!locNode->right) {
		avlRemoveFromParent(tree,locNode,locNode->left);
		locNode->left->parent = locNode->parent;
		if (locNode->parent)
		    avlUpdateMaxScores(locNode->parent);
		locNode->left = NULL;
		if (freeNodeMem)
		    avlFreeNode(ctx, locNode,0);
		*removed = 1;
		return -1;
	    }

	    // If two children, replace from subtree
	    if (locNode->balance < 0) {
		// Replace with the node's in-order predecessor
		replacementNode = locNode->left;
		while (replacementNode->right)
		    replacementNode = replacementNode->right;
	    }
	    else {
		// Replace with the node's in-order successor
		replacementNode = locNode->right;
		while (replacementNode->left)
		    replacementNode = replacementNode->left;
	    }

	    // Remove the replacementNode from the tree
	    heightDelta = avlRemoveNode(ctx, tree, locNode,replacementNode,0,removed);

	    //if the replacement node is a direct child of locNode,
	    //don't set it up to point to itself
	    if (locNode->right) locNode->right->parent = replacementNode;
	    if (locNode->left)  locNode->left->parent = replacementNode;

	    replacementNode->left = locNode->left;
	    replacementNode->right = locNode->right;
	    replacementNode->parent = locNode->parent;
	    replacementNode->balance = locNode->balance;

	    //Now replace the tree's root with replacementNode if it's the root
	    //otherwise place the replacement under the parent
	    if (locNode == tree->root) {
		tree->root = replacementNode;

		avlUpdateMaxScores(replacementNode);
	    }
	    else {
		if (locNode == locNode->parent->left)
		    locNode->parent->left = replacementNode;
		else
		    locNode->parent->right = replacementNode;

		avlUpdateMaxScores(replacementNode);
	    }

	    locNode->left = NULL;
	    locNode->right = NULL;
	    if (freeNodeMem)
		avlFreeNode(ctx, locNode,0);

	    if (replacementNode->balance == 0)
		return heightDelta;

	    *removed = 1;

	    return 0;
	}
    }

    // The node is in the left subtree
    else if (diff > 0) {
	if (locNode->left) {
	    heightDelta = avlRemoveNode(ctx,tree,locNode->left,delNode,freeNodeMem,removed);
	    if (heightDelta) {
		locNode->balance = locNode->balance + 1;
		if (locNode->balance == 0)
		    return -1;
		else if (locNode->balance == 1)
		    return 0;

		if (locNode->right->balance == 1) {
		    avlLeftRotation(tree,locNode);
		    locNode->parent->balance = 0;
		    locNode->parent->left->balance = 0;

		    locNode->subRightMax = locNode->parent->subLeftMax;
		    locNode->parent->subLeftMax = -INFINITY;
		    avlUpdateMaxScores(locNode->parent);

		    return -1;
		}
		else if (locNode->right->balance == 0){
		    avlLeftRotation(tree,locNode);
		    locNode->parent->balance = -1;
		    locNode->parent->left->balance = 1;

		    locNode->subRightMax = locNode->parent->subLeftMax;
		    locNode->parent->subLeftMax = -INFINITY;
		    avlUpdateMaxScores(locNode->parent);

		    return 0;
		}
		avlRightRotation(tree,locNode->right);
		avlLeftRotation(tree,locNode);
		avlResetBalance(locNode->parent);

		locNode->subRightMax = locNode->parent->subLeftMax;
		locNode->parent->right->subLeftMax = locNode->parent->subRightMax;
		locNode->parent->subRightMax = -INFINITY;
		locNode->parent->subLeftMax = -INFINITY;

		avlUpdateMaxScores(locNode->parent);

		return -1;
	    }
	}
    }

    // The node is in the right subtree
    else if (diff < 0) {
	if (locNode->right) {
	    heightDelta = avlRemoveNode(ctx,tree,locNode->right,delNode,freeNodeMem,removed);
	    if (heightDelta) {
		locNode->balance = locNode->balance - 1;
		if (locNode->balance == 0)
		    return 1;
		else if (locNode->balance == -1)
		    return 0;

		if (locNode->left->balance == -1) {
		    avlRightRotation(tree,locNode);
		    locNode->parent->balance = 0;
		    locNode->parent->right->balance = 0;

		    locNode->subLeftMax = locNode->parent->subRightMax;
		    locNode->parent->subRightMax = -INFINITY;
		    avlUpdateMaxScores(locNode->parent);

		    return -1;
		}
		else if (locNode->left->balance == 0){
		    avlRightRotation(tree,locNode);
		    locNode->parent->balance = 1;
		    locNode->parent->right->balance = -1;

		    locNode->subLeftMax = locNode->parent->subRightMax;
		    locNode->parent->subRightMax = -INFINITY;
		    avlUpdateMaxScores(locNode->parent);

		    return 0;
		}
		avlLeftRotation(tree,locNode->left);
		avlRightRotation(tree,locNode);
		avlResetBalance(locNode->parent);

		locNode->subLeftMax = locNode->parent->subRightMax;
		locNode->parent->left->subRightMax = locNode->parent->subLeftMax;

		locNode->parent->subRightMax = -INFINITY;
		locNode->parent->subLeftMax = -INFINITY;
		avlUpdateMaxScores(locNode->parent);

		return -1;
	    }
	}
    }

    return 0;
}

int avlRemove(RedisModuleCtx *ctx, avl *tree, double lscore, double rscore, RedisModuleString * obj) {
    int removed = 0;

    if (!tree->root)
	return 0;

    avlNode *delNode = avlCreateNode(ctx, lscore, rscore, obj);
    avlRemoveNode(ctx, tree, tree->root, delNode, 1, &removed);
    avlFreeNode(ctx,delNode,0);

    if (removed)
	tree->size = tree->size - 1;

    if (tree->size == 0)
	tree->root = NULL;

    return removed;
}

size_t isetLength(const void *obj) {
    return ((avl *) obj)->size;
}

/*
This structure is a simple linked list that is built during the
avlFind method. On each successful stab, the stabbed node is
added to the list and the list tail is updated to the new node.
The genericStabCommand maintains a pointer to the list head, which
is the initial node passed to avlFind.
*/
typedef struct avlResultNode {
    avlNode * data;
    struct avlResultNode * next;
} avlResultNode;

avlResultNode *avlCreateResultNode(avlNode *data) {
    avlResultNode *arn = RedisModule_Alloc(sizeof(*arn));
    arn->data = data;
    arn->next = NULL;
    return arn;
}

void avlFreeResults(avlResultNode *node) {
    if (node->next)
	avlFreeResults(node->next);
    RedisModule_Free(node);
}

avlResultNode * avlStab(avlNode *node, double min, double max, avlResultNode *results) {

    // If the minimum endpoint of the interval falls to the right of the current node's interval and
    // any sub-tree intervals, there cannot be an interval match
    if (min > node->subRightMax && min > node->subLeftMax && min > node->scores[1])
	return results;

    // Search the node's left subtree
    if (node->left)
	results = avlStab(node->left, min, max, results);

    // Check to see if this node overlaps.
    // For now we're only going to check for containment.
    if (min >= node->scores[0] && max <= node->scores[1]) {
	avlResultNode * newResult = avlCreateResultNode(node);
	newResult->next = results;
	results = newResult;
    }

    // If the max endpoint of the interval falls to the left of the start of the current node's
    // interval, there cannot be an interval match to the right of this node
    if (max < node->scores[0])
	return results;

    // Search the node's right subtree
    if (node->right)
	results = avlStab(node->right, min, max, results);

    return results;
}

/*-----------------------------------------------------------------------------
 * Interval set commands
 *----------------------------------------------------------------------------*/



int iaddCommand(RedisModuleCtx *ctx, RedisModuleString **argv, int argc) {
    RedisModule_AutoMemory(ctx);
    RedisModuleString *key = argv[1];
    RedisModuleString *ele;
    double * curscores;
    double min = 0, max = 0;
    int added = 0;
    double *mins, *maxes;
    int j, elements = (argc-2)/3;
    avl *tree;
    avlNode *addedNode;
    dictEntry *de;
    /* 5, 8, 11... arguments */
    if ((argc - 2) % 3) {
	RedisModule_WrongArity(ctx);
	return REDISMODULE_ERR;
    }

    /* Start parsing all the scores, we need to emit any syntax error
     * before executing additions to the sorted set, as the command should
     * either execute fully or nothing at all. */
    mins  = RedisModule_Alloc(sizeof(double)*elements);
    maxes = RedisModule_Alloc(sizeof(double)*elements);
    for (j = 0; j < elements; j++) {
	/* mins are 2, 5, 8... */
	if (RedisModule_StringToDouble(argv[2+j*3],&mins[j])
	    != REDISMODULE_OK)
	{
	    RedisModule_Free(mins);
	    RedisModule_Free(maxes);
	    RedisModule_ReplyWithError(ctx,REDISMODULE_ERRORMSG_SYNTAXERR);
	    return REDISMODULE_ERR;
	}
	/* maxes are 3, 6, 9... */
	if (RedisModule_StringToDouble(argv[3+j*3],&maxes[j])
	    != REDISMODULE_OK)
	{
	    RedisModule_Free(mins);
	    RedisModule_Free(maxes);
	    RedisModule_ReplyWithError(ctx,REDISMODULE_ERRORMSG_SYNTAXERR);
	    return REDISMODULE_ERR;
	}
    }

    /* Lookup the key and create the interval tree if does not exist. */
    RedisModuleKey *ikey = RedisModule_OpenKey(ctx,key,REDISMODULE_READ|REDISMODULE_WRITE);
    if (RedisModule_KeyType(ikey) != REDISMODULE_KEYTYPE_EMPTY &&
	RedisModule_ModuleTypeGetType(ikey) != ISet)
      {
	    RedisModule_Free(mins);
	    RedisModule_Free(maxes);
	    RedisModule_ReplyWithError(ctx,REDISMODULE_ERRORMSG_WRONGTYPE);
	    return REDISMODULE_ERR;
    }

    tree = RedisModule_ModuleTypeGetValue(ikey);

    if (tree == NULL) {
      tree = avlCreate();
      if(RedisModule_ModuleTypeSetValue(ikey,ISet,tree) != REDISMODULE_OK) {
	  return REDISMODULE_ERR;
	}
    }
    for (j = 0; j < elements; j++) {
	min = mins[j];
	max = maxes[j];

	ele = argv[4+j*3];

	de  = dictFind(tree->dict,ele);
	/* If object is found in iobj */

	if (de != NULL) {
	    curscores = (double *) dictGetVal(de);

	    if (curscores[0] != min || curscores[1] != max) {
		// remove and re-insert
		avlRemove(ctx, tree, curscores[0], curscores[1], ele);
		addedNode = avlInsert(ctx, tree, min, max, ele);
		dictGetVal(de) = addedNode->scores; /* Update scores ptr. */
		added++;
	    }
	} else {
	    /* insert into the tree */
	    addedNode = avlInsert(ctx, tree, min, max, ele);
	    if(!(dictAdd(tree->dict,ele,&addedNode->scores) == DICT_OK)) return REDISMODULE_ERR;
	    RedisModule_RetainString(ctx,ele); /* Added to dictionary. */
	    added++;
	}
    }

    RedisModule_Free(mins);
    RedisModule_Free(maxes);

    RedisModule_CloseKey(ikey);
    RedisModule_ReplyWithLongLong(ctx,added);

    return REDISMODULE_OK;
}

/* This command implements ISTAB, ISTABINTERVAL. */
int genericStabCommand(RedisModuleCtx *ctx, RedisModuleString *lscoreObj, RedisModuleString *rscoreObj, int intervalstab, RedisModuleString **argv, int argc) {
    RedisModule_AutoMemory(ctx);
    double lscore, rscore;
    int withintervals = 0;
    RedisModuleString *key = argv[1];
    avlResultNode * resnode;
    avlResultNode * reswalker;
    avlNode * nodewalker;
    unsigned long resultslen = 0;
    avl * tree;

    if (intervalstab) {
	if (argc > 4) {
	    if (!RedisModule_StringCompare(argv[4],RedisModule_CreateString(ctx,"withintervals",14)))
		withintervals = 1;
	    else {
		RedisModule_ReplyWithError(ctx,REDISMODULE_ERRORMSG_SYNTAXERR);
		return REDISMODULE_ERR;
	    }
	}
    } else {
	if (argc > 3) {
	    if (!RedisModule_StringCompare(argv[3],RedisModule_CreateString(ctx,"withintervals",14)))
		withintervals = 1;
	    else {
		RedisModule_ReplyWithError(ctx,REDISMODULE_ERRORMSG_SYNTAXERR);
		return REDISMODULE_ERR;
	    }
	}
    }

    if (RedisModule_StringToDouble(lscoreObj,&lscore) != REDISMODULE_OK) {
	return REDISMODULE_ERR;
    }

    if (RedisModule_StringToDouble(rscoreObj,&rscore) != REDISMODULE_OK) {
	return REDISMODULE_ERR;
    }

    RedisModuleKey *ikey = RedisModule_OpenKey(ctx,key,REDISMODULE_READ|REDISMODULE_WRITE);
    if (RedisModule_KeyType(ikey) == REDISMODULE_KEYTYPE_EMPTY &&
	RedisModule_ModuleTypeGetType(ikey) != ISet)
      {
	    RedisModule_ReplyWithError(ctx,REDISMODULE_ERRORMSG_WRONGTYPE);
	    return REDISMODULE_ERR;
    }

    tree = RedisModule_ModuleTypeGetValue(ikey);

    resnode = avlStab(tree->root, lscore, rscore, NULL);

    /* No results. */
    if (resnode == NULL) {
	RedisModule_ReplyWithNull(ctx);
	return REDISMODULE_ERR;
    }

    /* We don't know in advance how many matching elements there are in the
     * list, so we push this object that will represent the multi-bulk
     * length in the output buffer, and will "fix" it later */
    RedisModule_ReplyWithArray(ctx,REDISMODULE_POSTPONED_ARRAY_LEN);
    reswalker = resnode;

    while (reswalker != NULL) {
	nodewalker = reswalker->data;
	while (nodewalker != NULL) {
	    resultslen++;
	    RedisModule_ReplyWithString(ctx,nodewalker->obj);
	    /* if (withintervals) { */
	    /*	RedisModule_ReplyWithDouble(ctx,nodewalker->scores[0]); */
	    /*	RedisModule_ReplyWithDouble(ctx,nodewalker->scores[1]); */
	    /* } */
	    nodewalker = nodewalker->next;
	}
	reswalker = reswalker->next;
    }

    if (withintervals)
	resultslen *= 3;

    RedisModule_ReplySetArrayLength(ctx,resultslen);

    if (resnode)
	avlFreeResults(resnode);
    return REDISMODULE_OK;
}

int istabCommand(RedisModuleCtx *ctx,RedisModuleString **argv, int argc) {
    return genericStabCommand(ctx,argv[2],argv[2],0,argv,argc);
}

int istabIntervalCommand(RedisModuleCtx *ctx,RedisModuleString **argv, int argc) {
    return genericStabCommand(ctx,argv[2],argv[3],1,argv,argc);
}

void irembystabCommand(RedisModuleCtx *ctx,RedisModuleString **argv) {
    RedisModuleString *key = argv[1];
    long long point;
    avlResultNode * resnode;
    avlResultNode * reswalker;
    avlNode * nodewalker;
    avl * tree;
    int deleted = 0;
    dictEntry *de;
    double *curscores;

    if (RedisModule_StringToLongLong(argv[2], &point) != REDISMODULE_OK) return;

    /* Lookup the key and create the interval tree if does not exist. */
    RedisModuleKey *ikey = RedisModule_OpenKey(ctx,key,REDISMODULE_READ|REDISMODULE_WRITE);
    int type = RedisModule_KeyType(ikey);
    if (type != REDISMODULE_KEYTYPE_EMPTY &&
	RedisModule_ModuleTypeGetType(ikey) != ISet)
      {
	    RedisModule_ReplyWithError(ctx,REDISMODULE_ERRORMSG_WRONGTYPE);
	    return;
    }

    tree = RedisModule_ModuleTypeGetValue(ikey);

    resnode = avlStab(tree->root, point, point, NULL);

    /* No results. */
    if (resnode == NULL) {
	RedisModule_ReplyWithLongLong(ctx, 0);
	return;
    }

    reswalker = resnode;

    while (reswalker != NULL) {
	nodewalker = reswalker->data;
	while (nodewalker != NULL) {
	    de = dictFind(tree->dict,nodewalker->obj);
	    if (de != NULL) {
		deleted++;

		/* delete from the tree */
		curscores = (double *) dictGetVal(de);
		avlRemove(ctx,tree,curscores[0],curscores[1],nodewalker->obj);

		/* delete from the hash table */
		dictDelete(tree->dict,nodewalker->obj);

		if (dictSize(tree->dict) == 0) {
		    RedisModule_UnlinkKey(ikey);
		    break;
		}
	    }
	    nodewalker = nodewalker->next;
	}
	reswalker = reswalker->next;
    }

    if (resnode)
	avlFreeResults(resnode);

    RedisModule_ReplyWithLongLong(ctx,deleted);
}

int iremCommand(RedisModuleCtx *ctx,RedisModuleString **argv, int argc) {
    RedisModule_AutoMemory(ctx);
    RedisModuleString *key = argv[1];
    RedisModuleString *ele;
    double *curscores;
    int deleted = 0, j;
    dictEntry *de;
    avl *tree;

    /* Lookup the key and create the interval tree if does not exist. */
    RedisModuleKey *ikey = RedisModule_OpenKey(ctx,key,REDISMODULE_READ|REDISMODULE_WRITE);
    int type = RedisModule_KeyType(ikey);
    if (type != REDISMODULE_KEYTYPE_EMPTY &&
	RedisModule_ModuleTypeGetType(ikey) != ISet)
      {
	    RedisModule_ReplyWithError(ctx,REDISMODULE_ERRORMSG_WRONGTYPE);
	    return REDISMODULE_ERR;
    }

    tree = RedisModule_ModuleTypeGetValue(ikey);

    for (j = 2; j < argc; j++) {
	ele = argv[j];
	de = dictFind(tree->dict,ele);

	if (de != NULL) {
	    deleted++;

	    /* delete from the tree */
	    curscores = (double *) dictGetVal(de);
	    avlRemove(ctx,tree,curscores[0],curscores[1],ele);

	    /* delete from the hash table */
	    dictDelete(tree->dict,ele);
	    if (dictSize(tree->dict) == 0) {
		RedisModule_UnlinkKey(ikey);
		break;
	    }
	}
    }

    RedisModule_ReplyWithLongLong(ctx,deleted);
    return REDISMODULE_OK;
}

void *isetRdbLoad(RedisModuleIO *rdb, int encver){
  REDISMODULE_NOT_USED(encver);
  uint64_t isetlen;
  avl *tree;
  RedisModuleCtx *ctx;

  ctx     = RedisModule_GetContextFromIO(rdb);
  tree    = avlCreate();
  isetlen = RedisModule_LoadUnsigned(rdb);

  while(isetlen--) {
    RedisModuleString *ele;
    double score1;
    double score2;

    avlNode *inode;

    ele    = RedisModule_LoadString(rdb);
    score1 = RedisModule_LoadDouble(rdb);
    score2 = RedisModule_LoadDouble(rdb);

    inode = avlInsert(ctx, tree, score1, score2, ele);
    dictAdd(tree->dict,ele,&inode->scores);
    RedisModule_RetainString(ctx,ele); /* Added to dictionary. */
  }
  return avlCreate;
}

void isetRdbSave(RedisModuleIO *rdb, void *value){
  avl * tree = value;
  dictIterator *di = dictGetIterator(tree->dict);
  dictEntry *de;
  uint64_t isetlen = isetLength(value);

  RedisModule_SaveUnsigned(rdb,isetlen);
  while((de = dictNext(di)) != NULL) {
    RedisModuleString *eleobj = dictGetKey(de);
    double *scores = dictGetVal(de);
    RedisModule_SaveString(rdb,eleobj);
    RedisModule_SaveDouble(rdb,scores[0]);
    RedisModule_SaveDouble(rdb,scores[1]);
  }
  dictReleaseIterator(di);
}

void rewriteIntervalSetObject(RedisModuleIO *r, RedisModuleString *key, void *o) {
    dictIterator *di = dictGetIterator(((avl *) o)->dict);
    dictEntry *de;
    const char *keyStr = RedisModule_StringPtrLen(key,NULL);

    while((de = dictNext(di)) != NULL) {
      RedisModuleString *eleobj = dictGetKey(de);
      double **scores           = dictGetVal(de);
      const char *eleStr        = RedisModule_StringPtrLen(eleobj,NULL);
      RedisModule_EmitAOF(r,"IADD","%s %d %d %s", keyStr, (*scores)[0], (*scores)[1], eleStr);
    }
    dictReleaseIterator(di);
}

int RedisModule_OnLoad(RedisModuleCtx *ctx, RedisModuleString **argv, int argc) {
    REDISMODULE_NOT_USED(argv);
    REDISMODULE_NOT_USED(argc);

    if (RedisModule_Init(ctx,"intervalS",1,REDISMODULE_APIVER_1)
	== REDISMODULE_ERR) return REDISMODULE_ERR;

    RedisModuleTypeMethods tm = {
	.version     = REDISMODULE_TYPE_METHOD_VERSION,
	/* .rdb_load    = isetRdbLoad, */
	/* .rdb_save    = isetRdbSave, */
	/* .aof_rewrite = rewriteIntervalSetObject, */
	.free        = avlFree,
	/* .mem_usage   = isetLength, */
	/* .digest      = NULL */
    };

    ISet = RedisModule_CreateDataType(ctx,"intervalS",0,&tm);
    if (ISet == NULL) return REDISMODULE_ERR;
    if (RedisModule_CreateCommand(ctx,"intervallSet.insert",
	iaddCommand,"write deny-oom",1,1,1) == REDISMODULE_ERR)
	return REDISMODULE_ERR;

    if (RedisModule_CreateCommand(ctx,"intervallSet.rem",
	iremCommand,"write",1,1,1) == REDISMODULE_ERR)
	return REDISMODULE_ERR;

    if (RedisModule_CreateCommand(ctx,"intervallSet.stab",
	istabCommand,"write",1,1,1) == REDISMODULE_ERR)
	return REDISMODULE_ERR;

    if (RedisModule_CreateCommand(ctx,"intervallSet.stabIntervall",
	istabIntervalCommand,"write",1,1,1) == REDISMODULE_ERR)
	return REDISMODULE_ERR;


    return REDISMODULE_OK;
}
