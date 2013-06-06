#ifndef __NEWICKFORM_H__
#define __NEWICKFORM_H__
#include <stdbool.h>

typedef struct newick_child
{
	struct newick_node *node;
	struct newick_child *next;
} newick_child;

typedef struct newick_node
{
	char *taxon;
	char *seq;
	float dist;
	int childNum;
	struct newick_child *child;
	struct newick_node *parent;
} newick_node;

typedef struct sAllLeafs
{
	newick_child** 	paChildren;		/* Pointer to array, needs to be malloc'ed */
	int 			curIndex;  	/*  Current position */
	int 			size;			/*  Leaf count */
} sAllLeafs;

short sDEBUG;

#ifdef __NEWICKFORM_C__
newick_node* parseTree(char *str);
void printTree(newick_node *root);
void traverseStructure(newick_node *root, newick_child **newick);
void storeTaxa(newick_node *root);
void countTaxaNumber(newick_node *root, int* i);
bool leafs_add_leaf(sAllLeafs* all_leafs, newick_child* c);

#else
extern newick_node* parseTree(char *str);
extern void printTree(newick_node *root);
extern void traverseStructure(newick_node *root, newick_child **newick);
extern void storeTaxa(newick_node *root);
extern void countTaxaNumber(newick_node *root, int* i);
extern bool leafs_add_leaf(sAllLeafs* all_leafs, newick_child* c);

#endif

#endif

