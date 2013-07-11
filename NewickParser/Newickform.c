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
			printf(")%s:%0.5f", root->taxon, root->dist);
		}
		else
		{
			printf("):%0.5f", root->dist);
		}
	}
}


/* return leafs for given index */
newick_node* leafs_get_leaf(sAllLeafs* all_leafs, int index) {
	newick_node* result = NULL;

	if (all_leafs != NULL && index < all_leafs->size) {
		result = all_leafs->paChildren[index];
	}

	return result;
}

/* Add a leaf to a given sAllLeafs* pointer */
bool leafs_add_leaf(sAllLeafs* all_leafs, newick_node* c) {
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
	all_leafs->paChildren = (newick_node**)malloc(leaf_count * sizeof(newick_node*));
	all_leafs->curIndex = 0;
	all_leafs->size = leaf_count;

	return all_leafs;
}

void generateStringFromArray(sAllLeafs* arr, int i) {
    printf("%s:%0.5f", leafs_get_leaf(arr, i)->taxon, leafs_get_leaf(arr, i)->dist);
}

/* Store all Taxa in an array called all_leafs -- TESTFUNCTION()*/
void storeTaxa(newick_node *root){
    int TaxaNumber = 0;
    countTaxaNumber(root, &TaxaNumber);
    sAllLeafs* all_leafs = leafs_alloc(TaxaNumber);
    traverseStructure(root, &all_leafs);

    printf("Subtrees: \n");
    //int i = 0;
    //leafs_get_leaf(all_leafs, 2)->taxon = "human";
    randomize(all_leafs, TaxaNumber);
    //while (i < TaxaNumber) {
        
        //generateStringFromArray(all_leafs, i);
         //i++;
    //}
    parens(all_leafs, all_leafs->size / 4);
    //printf("Size: %d \n", all_leafs->size);
    
    
}

/* Count Taxa Number */
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
	}

}



void generate_randomTree(int leafs) {
    newick_node *root;
	newick_child *child;
    
    root = (newick_node*)seqMalloc(sizeof(newick_node));
    child = (newick_child*)seqMalloc(sizeof(newick_child));
    
    
    for (int i = 0; i < leafs; i++) {
        
    }
}

/*returns all taxa */
void traverseStructure(newick_node *root, sAllLeafs** newick)
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
        //printf("%s, %f \n", root->taxon, root->dist);
        leafs_add_leaf(*newick, root);
	}

}
//Swap nodes in position a and b (beginning with 0)
void leafs_swap_nodes(int a, int b, sAllLeafs* i)
{
    newick_node* tmp = leafs_get_leaf(i, a);
    i->paChildren[a] = leafs_get_leaf(i, b);
    i->paChildren[b] = tmp;
}

//Randomizes array
void randomize(sAllLeafs* arr, int n )
{
// Use a different seed value so that we don't get same
// result each time we run this program
srand ( time(NULL) );
    
// Start from the last element and swap one by one. We don't
// need to run for the first element that's why i > 0
    for (int i = n-1; i > 0; i--)
        {  
            // Pick a random index from 0 to i
            int j = rand() % (i+1);
            
            // Swap arr[i] with the element at random index
            leafs_swap_nodes(i, j, arr);
        }
}


void parens_foo(char output[], int index, int open, int close, int pairs, sAllLeafs* arr)
{
    int i;
    int k = 0;
    
    if (index == 2*pairs) {
        printf("(");
        for (i = 0; i < 2*pairs; i++) {
            putchar(output[i]);
            if (output[i] == '(') {
                generateStringFromArray(arr, k);
                k++;
                if(output[i+1] == ')') {
                    printf(",");
                    generateStringFromArray(arr, k);
                    k++;
                }
                if(output[i+1] == '(') {
                    printf(",");
                }
            }
            if (output[i] == ')') {
                printf(":10.00000");
                if (output[i+1] == '(') {
                    printf(",");
                }
            }
            
        }
        printf(",");
        generateStringFromArray(arr, k);
        printf(",");
        k++;
        generateStringFromArray(arr, k);
        printf(");");
        putchar('\n');
        
        return;
    }
    
    if (open) {
        output[index] = '(';
        parens_foo(output, index + 1, open - 1, close, pairs, arr);
    }
    
    if (close && (pairs - close + 1 <= pairs - open)
        ) {
        output[index] = ')';
        parens_foo(output, index + 1, open, close - 1, pairs, arr);
    }
    
    return;
}

//Generate Parenthesis
void parens(sAllLeafs* arr, int pairs)
{
    char output[2*pairs];
    
    output[0] = '(';
     output[2*pairs - 1] = ')';
    
    parens_foo(output, 1, pairs - 1, pairs, pairs, arr);
}








