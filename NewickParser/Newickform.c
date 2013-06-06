#define __NEWICKFORM_C__

#include "seqUtil.h"
#include "Newickform.h"

newick_node* parseTree(char *str)
{
	newick_node *node;
	newick_child *child;
	char *pcCurrent;
	char *pcStart;
	char *pcColon = NULL;
	char cTemp;
	int iCount;

	if (sDEBUG == 1)
	{
		printf("%s\n", str);
	}
	pcStart = str;

	if (*pcStart != '(')
	{
		// Leaf node. Separate taxon name from distance. If distance not exist then take care of taxon name only
		pcCurrent = str;
		while (*pcCurrent != '\0')
		{
			if (*pcCurrent == ':')
			{
				pcColon = pcCurrent;
			}
			pcCurrent++;
		}
		node = (newick_node*)seqMalloc(sizeof(newick_node));
		if (pcColon == NULL)
		{
			// Taxon only
			node->taxon = (char*)seqMalloc(strlen(pcStart) + 1);
			memcpy(node->taxon, pcStart, strlen(pcStart));
		}
		else
		{
			// Taxon
			*pcColon = '\0';
			node->taxon = (char*)seqMalloc(strlen(pcStart) + 1);
			memcpy(node->taxon, pcStart, strlen(pcStart));
			*pcColon = ':';
			// Distance
			pcColon++;
			node->dist = (float)atof(pcColon);
		}
		node->childNum = 0;
	
    }
	else
	{
		// Create node
		node = (newick_node*)seqMalloc(sizeof(newick_node));
		child = NULL;
		// Search for all child nodes
		// Find all ',' until corresponding ')' is encountered
		iCount = 0;
		pcStart++;
		pcCurrent = pcStart;
		while (iCount >= 0)
		{
			switch (*pcCurrent)
			{
				case '(':
					// Find corresponding ')' by counting
					pcStart = pcCurrent;
					pcCurrent++;
					iCount++;
					while (iCount > 0)
					{
						if (*pcCurrent == '(')
						{
							iCount++;
						}
						else if (*pcCurrent == ')')
						{
							iCount--;
						}
						pcCurrent++;
					}
					while (*pcCurrent != ',' && *pcCurrent != ')')
					{
						pcCurrent++;
					}
 					cTemp = *pcCurrent;
					*pcCurrent = '\0';
					// Create a child node
					if (child == NULL)
					{
						node->child = (newick_child*)seqMalloc(sizeof(newick_child));
						node->childNum = 1;
						child = node->child;
					}
					else
					{
						child->next = (newick_child*)seqMalloc(sizeof(newick_child));
						node->childNum++;
						child = child->next;
					}
					child->node = parseTree(pcStart);
					*pcCurrent = cTemp;
					if (*pcCurrent != ')')
					{
						pcCurrent++;
					}
				break;

				case ')':
					// End of tihs tree. Go to next part to retrieve distance
					iCount--;
				break;

				case ',':
					// Impossible separation since according to the algorithm, this symbol will never encountered.
					// Currently don't handle this and don't create any node
				break;

				default:
					// leaf node encountered
					pcStart = pcCurrent;
					while (*pcCurrent != ',' && *pcCurrent != ')')
					{
						pcCurrent++;
					}
					cTemp = *pcCurrent;
					*pcCurrent = '\0';
					// Create a child node
					if (child == NULL)
					{
						node->child = (newick_child*)seqMalloc(sizeof(newick_child));
						node->childNum = 1;
						child = node->child;
					}
					else
					{
						child->next = (newick_child*)seqMalloc(sizeof(newick_child));
						node->childNum++;
						child = child->next;
					}
					child->node = parseTree(pcStart);
					*pcCurrent = cTemp;
					if (*pcCurrent != ')')
					{
						pcCurrent++;
					}
				break;
			}
		}

		// If start at ':', then the internal node has no name.
		pcCurrent++;
		if (*pcCurrent == ':')
		{
			pcStart = pcCurrent + 1;
			while (*pcCurrent != '\0' && *pcCurrent != ';')
			{
				pcCurrent++;
			}
			cTemp = *pcCurrent;
			*pcCurrent = '\0';
			node->dist = (float)atof(pcStart);
			*pcCurrent = cTemp;
		}
		else if (*pcCurrent != ';' && *pcCurrent != '\0')
		{
			// Find ':' to retrieve distance, if any.
			// At this time *pcCurrent should equal to ')'
			pcStart = pcCurrent;
			while (*pcCurrent != ':')
			{
				pcCurrent++;
			}
			cTemp = *pcCurrent;
			*pcCurrent = '\0';
			node->taxon = seqMalloc(strlen(pcStart) + 1);
			memcpy(node->taxon, pcStart, strlen(pcStart));
			*pcCurrent = cTemp;
			pcCurrent++;
			pcStart = pcCurrent;
			while (*pcCurrent != '\0' && *pcCurrent != ';')
			{
				pcCurrent++;
			}
			cTemp = *pcCurrent;
			*pcCurrent = '\0';
			node->dist = (float)atof(pcStart);
			*pcCurrent = cTemp;
		}
	}

	return node;
}

void printTree(newick_node *root)
{
	newick_child *child;
	/* Check if node is a tip */
	if (root->childNum == 0)
	{
		printf("%s:%0.6f", root->taxon, root->dist);
	}
	else
	{
		child = root->child;
		printf("(");
		while (child != NULL)
		{
			printTree(child->node);
			if (child->next != NULL)
			{
				printf(",");
			}
			child = child->next;
		}
		if (root->taxon != NULL)
		{
			printf(")%s:%0.6f", root->taxon, root->dist);
		}
		else
		{
			printf("):%0.6f", root->dist);
		}
	}
}



newick_child* leafs_get_leaf(sAllLeafs* all_leafs, int index) {
	newick_child* result = NULL;

	if (all_leafs != NULL && index < all_leafs->size) {
		result = all_leafs->paChildren[index];
	}

	return result;
}

/* Add a leaf to a given sAllLeafs* pointer */
bool leafs_add_leaf(sAllLeafs* all_leafs, newick_child* c) {
	bool result = false;

	/* Check that leaf structure is not null */
	if (all_leafs != NULL && all_leafs->curIndex < all_leafs->size) {
		/* Add at current position */
		all_leafs->paChildren[all_leafs->curIndex] = c;
		/* Increment current position */
		all_leafs->curIndex++;
		result = true;
	}

	return result;
}

/* Destructor - freeing memory */
void leafs_free(sAllLeafs* all_leafs) {
	free(all_leafs);
}


/* Constructor - allocating memory, setting struct values */
sAllLeafs* leafs_alloc(int leaf_count) {
	/* Example of allocation */
	sAllLeafs* all_leafs = (sAllLeafs*)malloc(sizeof(sAllLeafs));
	all_leafs->paChildren = (newick_child**)malloc(leaf_count * sizeof(newick_child*));
	all_leafs->curIndex = 0;
	all_leafs->size = leaf_count;

	return all_leafs;
}




void storeTaxa(newick_node *root){
int TaxaNumber = 0;
newick_child* child;
child = (newick_child*)seqMalloc(sizeof(newick_child));
countTaxaNumber(root, &TaxaNumber);
int i = 0;
while(i <= TaxaNumber) {
	child->next = (newick_child*)seqMalloc(sizeof(newick_child));
	child = child->next;
	i++;
}
traverseStructure(root, &child); 

printf("hey\n");
while(child != NULL) {
	printf("%s ", child->node->taxon);
	child = child->next;
	}
}

void countTaxaNumber(newick_node *root, int* i) {
	if(root->childNum != 0) {
		newick_child *child;
		child = root->child;
		if(child != NULL) {
		countTaxaNumber(child->node, i);
			if(child->next != NULL) {
			countTaxaNumber(child->next->node, i);
				if(child->next->next != NULL) {
				countTaxaNumber(child->next->next->node, i);
				}
			}

		}
	} else {
	(*i) = (*i) + 1;
	//printf("%i", (*i));
	//printf("%s \n", root->taxon);
	}

}


/*returns all taxa */
void traverseStructure(newick_node *root, newick_child** newick) 
{	
	if(root->childNum != 0) {
		newick_child *child;
		child = root->child;
		if(child != NULL) {
			traverseStructure(child->node, newick);
			if(child->next != NULL) {
				traverseStructure(child->next->node, newick);
				if(child->next->next != NULL) {
					traverseStructure(child->next->next->node, newick);
				}
			}	
		}
	} else {
	printf("%s \n", root->taxon);

	(*newick)->node = root; //NOT WORKING RIGHT NOW 
	newick = &((*newick)->next);
	}

}



