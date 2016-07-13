// Life.cpp : Defines the entry point for the console application.
//
//	This is an adaptation of HashLife Java implementation by Tomas G. Rokicki:
//	http://www.ddj.com/184406478
//  Modified to allow producing arbitrary number of steps, converted to C++
//
//	Name: Aliaksei Sanko
//	Blog: aliakseis.livejournal.com
//	Country: Belarus

#include "stdafx.h"

#include <iostream>
#include <fstream>

#include <shlwapi.h>
#pragma comment(lib, "shlwapi")

#include <assert.h>


using std::cerr;
using std::ofstream;
using std::endl;


enum { HASH_SIZE = 64 * 1024 };

/**
*   This class contains the tree maintenance functions for quadtrees.
*/
class Node 
{
public:
	/**
	*   Construct a leaf cell.
	*/
	Node(int living) 
	{
		nw = ne = sw = se = 0 ;
		level = 1;//0 ;
		alive = living ;

		result[0] = 0;
		result[1] = 0;
	}
	/**
	*   Construct a node given four children.
	*/
	Node(Node* nw_, Node* ne_, Node* sw_, Node* se_, int living) 
	{
		nw = nw_ ;
		ne = ne_ ;
		sw = sw_ ;
		se = se_ ;
		level = nw_->level + 1 ;
		alive = living;

		result[0] = 0;
		result[1] = 0;
	}
	/**
	*   Set a bit in this node in its relative coordinate system;
	*   returns a whole new node since our nodes are immutable.
	*
	*   In the recursive call, we simply adjust the coordinate system
	*   and call down a level.
	*/
	Node* setBit(int x, int y) {
		//if (level == 0)
		//	return &aliveNode;
		if (level == 1)
		{
			int idx = this - level1Nodes;
			idx |= 1 << (-y * 2 - x);
			return &level1Nodes[idx];
		}

		// distance from center of this node to center of subnode is
		// one fourth the size of this node.
		int offset = 1 << (level - 2) ;
		if (x < 0)
			if (y < 0)
				return create(nw->setBit(x+offset, y+offset), ne, sw, se) ;
			else
				return create(nw, ne, sw->setBit(x+offset, y-offset), se) ;
		else
			if (y < 0)
				return create(nw, ne->setBit(x-offset, y+offset), sw, se) ;
			else
				return create(nw, ne, sw, se->setBit(x-offset, y-offset)) ;
	}
	/**
	*   If we ever really need to get a bit one at a time, we can
	*   use this subroutine.  For convenience it returns 0/1 rather
	*   than false/true.
	*/
	int getBit(int x, int y) {

        if (!alive)
            return 0;
		//if (level == 0)
		//	return 1;
		if (level == 1)
			return (alive & (1 << (-y * 4 - x))) ? 1 : 0;
		int offset = 1 << (level - 2) ;
		if (x < 0)
			if (y < 0)
				return nw->getBit(x+offset, y+offset) ;
			else
				return sw->getBit(x+offset, y-offset) ;
		else
			if (y < 0)
				return ne->getBit(x-offset, y+offset) ;
			else
				return se->getBit(x-offset, y-offset) ;
	}
	/**
	*   Build an empty tree at the given level.
	*/
	static Node* emptyTree(int lev) {
		//if (lev == 0)
		//	return &emptyNode;
		if (lev == 1)
			return &level1Nodes[0];
		Node* n = emptyTree(lev-1) ;
		return create(n, n, n, n, false) ;
	}
	/**
	*   Expand the universe; return a new node up one level with the
	*   current node in the center.  Requires us to disassemble the
	*   current node.
	*/
	Node* expandUniverse() {
		Node* border = emptyTree(level-1) ;
		return create(create(border, border, border, nw, nw->alive),
			create(border, border, ne, border, ne->alive),
			create(border, sw, border, border, sw->alive),
			create(se, border, border, border, se->alive), alive) ;
	}

	/**
	*   HorizontalForward() takes two horizontally adjacent nodes,
	*   builds a new node from the east half of the west node and
	*   the west half of the east node, and computes the next
	*   step for that new node.
	*/
	Node* horizontalForward(Node* w, Node* e) {
		return create(w->ne, e->nw, w->se, e->sw)->nextGeneration() ;
	}
	/**
	*   VerticalForward() takes two vertically adjacent nodes,
	*   builds a new node from the south half of the north node
	*   and the north half of the south node, and computes the
	*   next step for that node.
	*/
	Node* verticalForward(Node* n, Node* s) {
		return create(n->sw, n->se, s->nw, s->ne)->nextGeneration() ;
	}
	/**
	*   CenterForward() builds a new subnode that is half the
	*   size of the this node out of the center portions.  It
	*   then does a generation step and returns that result.
	*/
	Node* centerForward() {
		return create(nw->se, ne->sw, sw->ne, se->nw)->nextGeneration() ;
	}


	/**
	*   Return a new node one level down containing only the
	*   center elements.
	*/
	Node* centeredSubnode() {
		return create(nw->se, ne->sw, sw->ne, se->nw) ;
	}
	/**
	*   Return a new node one level down from two given nodes
	*   that contains the east centered two sub sub nodes from
	*   the west node and the west centered two sub sub nodes
	*   from the east node.
	*/
	Node* centeredHorizontal(Node* w, Node* e) {
		return create(w->ne->se, e->nw->sw, w->se->ne, e->sw->nw) ;
	}
	/**
	*   Similar, but this does it north/south instead of east/west.
	*/
	Node* centeredVertical(Node* n, Node* s) {
		return create(n->sw->se, n->se->sw, s->nw->ne, s->ne->nw) ;
	}
	/**
	*   Return a new node two levels down containing only the
	*   centered elements.
	*/
	Node* centeredSubSubnode() {
		return create(nw->se->se, ne->sw->sw, sw->ne->ne, se->nw->nw) ;
	}

	/**
	*   NextGeneration() is the core HashLife recursive algorithm.
	*   It builds 9 subnodes that are one quarter the size of the
	*   current node and advanced in time one eighth the size of
	*   the current node.  It then takes these 9 subnodes in groups
	*   of four and builds 4 subnodes, each one quarter the size
	*   of the current node and advanced in time one fourth the size
	*   of the current node.  It combines these four subnodes into
	*   a new node and returns that as its result.
	*/
	Node* nextGeneration(unsigned long supplement = 0)
	{
		bool direct = level > 2 && (supplement & (1ul << (level - 3)));

		if (result[direct] != 0)
			return result[direct];
		if (!alive)
			return result[direct] = nw ;
		if (level == 2)
			return result[direct] = slowSimulation() ;

		Node *n00, *n01 ,*n02 ,*n10 ,*n11 ,*n12 ,*n20 ,*n21 ,*n22; 

		if (direct)
		{
			n00 = nw->centeredSubnode();
			n01 = centeredHorizontal(nw, ne);
			n02 = ne->centeredSubnode();
			n10 = centeredVertical(nw, sw);
			n11 = centeredSubSubnode();
			n12 = centeredVertical(ne, se);
			n20 = sw->centeredSubnode();
			n21 = centeredHorizontal(sw, se);
			n22 = se->centeredSubnode();
		}
		else
		{
			n00 = nw->nextGeneration();
			n01 = horizontalForward(nw, ne);
			n02 = ne->nextGeneration();
			n10 = verticalForward(nw, sw);
			n11 = centerForward();
			n12 = verticalForward(ne, se);
			n20 = sw->nextGeneration();
			n21 = horizontalForward(sw, se);
			n22 = se->nextGeneration();
		}

		return result[direct] = create(
			create(n00, n01, n10, n11)->nextGeneration(supplement),
			create(n01, n02, n11, n12)->nextGeneration(supplement),
			create(n10, n11, n20, n21)->nextGeneration(supplement),
			create(n11, n12, n21, n22)->nextGeneration(supplement)) ;
	}

	/**
	*   create functions.
	*/
	static Node* create(Node* nw, Node* ne, Node* sw, Node* se, int living = true) {
		return Node(nw, ne, sw, se, living).intern() ;
	}
	static Node* create() {
		return emptyTree(3) ;
	}

	/**
	*   We need to provide a hashCode() and an equals() method to be
	*   able to hash these objects.
	*/

	size_t hashCode() 
	{
        assert(level > 1);

		return (size_t(nw) +
			11 * size_t(ne) +
			101 * size_t(sw) +
			1007 * size_t(se)) >> 5;
	}
	
	bool operator ==(const Node& t) const
	{
        assert(level > 1);

		return nw == t.nw && ne == t.ne && sw == t.sw && se == t.se;
	}

	/**
	*   Given a node, return the canonical one if it exists, or make it
	*   the canonical one.
	*/
	Node* intern() 
	{
		int i =  hashCode() % HASH_SIZE;

		Node* canon;
		for (canon = hashTable[i]; canon != 0; canon = canon->next)
		{
			if (*this == *canon)
			{
				return canon;
			}
		}
		canon = new Node(*this);
		assert(canon->alive == (canon->ne->alive || canon->nw->alive || canon->se->alive || canon->sw->alive));
		canon->next = hashTable[i];
		hashTable[i] = canon;
		return canon;
	}

	/**
	*   Given an integer with a bitmask indicating which bits are
	*   set in the neighborhood, calculate whether this cell is
	*   alive or dead in the next generation.  The bottom three
	*   bits are the south neighbors; bits 4..6 are the current
	*   row with bit 5 being the cell itself, and bits 8..10
	*   are the north neighbors.
	*/
	//Node* oneGen(int bitmask) 
	//{
	//	int self = bitmask & (1 << 5);
	//	bitmask &= 0x757 ; // mask out bits we don't care about

	//	if (bitmask == 0)
	//		return &emptyNode;
	//	bitmask &= bitmask - 1 ; // clear least significant bit
	//	if (bitmask == 0)
	//		return &emptyNode;
	//	bitmask &= bitmask - 1 ; // clear least significant bit
	//	if (bitmask == 0)
	//		return self? &aliveNode : &emptyNode;
	//	bitmask &= bitmask - 1 ; // clear least significant bit
	//	return (bitmask == 0)? &aliveNode : &emptyNode;
	//}
	bool oneGen(int bitmask) 
	{
		int self = bitmask & (1 << 5);
		bitmask &= 0x757 ; // mask out bits we don't care about

		if (bitmask == 0)
			return 0;
		bitmask &= bitmask - 1 ; // clear least significant bit
		if (bitmask == 0)
			return 0;
		bitmask &= bitmask - 1 ; // clear least significant bit
		if (bitmask == 0)
			return self? 1 : 0;
		bitmask &= bitmask - 1 ; // clear least significant bit
		return (bitmask == 0)? 1 : 0;
	}


	/**
	*   At level 2, we can use slow simulation to compute the next
	*   generation.  We use bitmask tricks.
	*/
	int allBits()
	{
		if (alive)
			return (nw->alive << 5) + (ne->alive << 4) + (sw->alive << 1) + se->alive;
		return 0;
	}
	Node* slowSimulation() 
	{
		//int allbits = (nw->allBits() << 10) + (ne->allBits() << 8) + (sw->allBits() << 2) + se->allBits();
		int allbits = (nw->alive << 10) + (ne->alive << 8) + (sw->alive << 2) + se->alive;

		//return create(oneGen(allbits>>5), oneGen(allbits>>4),
		//	oneGen(allbits>>1), oneGen(allbits)) ;

		return &level1Nodes[(oneGen(allbits>>5) << 3) + (oneGen(allbits>>4) << 2) + (oneGen(allbits>>1) << 1) + oneGen(allbits)];
	}


    void* operator new(size_t count)
    {
        void* result = bufferPtr;
        bufferPtr += count;
        return result;
    }

    void operator delete(void*) 
    {
    }

	Node *nw, *ne, *sw, *se ; // our children
	int level ;           // distance to root
	int alive ;       // if leaf node, are we alive or dead?
	Node* result[2];
	Node* next;

//	static Node aliveNode;
//	static Node emptyNode;

	static Node level1Nodes[16];

    static char buffer[16 * 1024 * 1024]; 
    static char* bufferPtr;

	static Node* hashTable[HASH_SIZE];
};


//Node Node::aliveNode(true);
//Node Node::emptyNode(false);

Node Node::level1Nodes[16] = 
{
	0x00, 0x01, 0x02, 0x03,
	0x10, 0x11, 0x12, 0x13,
	0x20, 0x21, 0x22, 0x23,
	0x30, 0x31, 0x32, 0x33,
};

char Node::buffer[16 * 1024 * 1024]; 

char* Node::bufferPtr = Node::buffer;

Node* Node::hashTable[HASH_SIZE];


class Universe
{
public:
	Universe()
	{
		// Initialize static nodes stuff
		if (Node::bufferPtr != Node::buffer)
		{
			Node::bufferPtr = Node::buffer;
			memset(Node::hashTable, 0, sizeof(Node::hashTable));
		}

		generationCount = 0;
		root = Node::create();
	}

	void runSteps(unsigned long numSteps) 
	{
		//while (root->level < 3 ||
		//		root->nw->population != root->nw->se->se->population ||
		//		root->ne->population != root->ne->sw->sw->population ||
		//		root->sw->population != root->sw->ne->ne->population ||
		//		root->se->population != root->se->nw->nw->population)
		//	root = root->expandUniverse();

		unsigned long stepSize;
		while ((stepSize = 1ul << (root->level - 2)) < numSteps)
			root = root->expandUniverse();

		root = root->nextGeneration(stepSize - numSteps);
		generationCount += numSteps;
	}

	/**
	*   Set a single bit; can only do this before running, and once
	*   we've started running cannot change.
	*/
	void setBit(int x, int y) 
	{
		/*
		*   We need to make sure the current universe is large enough
		*   to handle these parameters.  A root node at level n supports
		*   coordinates from -2^(n-1) .. 2^(n-1)-1.
		*/
		for (;;) {
			int maxCoordinate = 1 << (root->level - 1) ;
			if (-maxCoordinate <= x && x <= maxCoordinate-1 &&
				-maxCoordinate <= y && y <= maxCoordinate-1)
				break ;
			root = root->expandUniverse() ;
		}
		/*
		*   Call our recursive routine to set the bit.
		*/
		root = root->setBit(x, y) ;
	}

	Node* root;
	unsigned long generationCount;
};

/////////////////////////////////////////////////////////////////////

void SetBit(Universe& universe, int x, int y)
{
	for (int i = -1000; i <= 1000; i += 1000)
		for (int j = -1000; j <= 1000; j += 1000)
			universe.setBit(x + i, y + j);
}


int main(int /*argc*/, char* argv[])
{
	LARGE_INTEGER frequency;
	QueryPerformanceFrequency(&frequency);
	LARGE_INTEGER start, stop;

	QueryPerformanceCounter(&start);
    
	int runsCount = 0;

	for (;;) 
	{
		Universe universe;

		SetBit(universe, 1, 0);
		SetBit(universe, 2, 0);
		SetBit(universe, 0, 1);
		SetBit(universe, 1, 1);
		SetBit(universe, 1, 2);

		universe.runSteps(1000);

		runsCount++;

		QueryPerformanceCounter(&stop);
		if ((stop.QuadPart - start.QuadPart) / frequency.QuadPart < 5)
			continue;

		char path[_MAX_PATH];
		strcpy(path, argv[0]);
		char* pFileName = PathFindFileNameA(path);

		strcpy(pFileName, "life.txt");

		ofstream outputFile(path);
		if (!outputFile) 
		{
			cerr << "Unable to open output file.\n";
			return EXIT_FAILURE;
		}

		for (int y = -500; y < 500; ++y)
		{
			for (int x = -500; x < 500; ++x)
			{
				outputFile << ((int) universe.root->getBit(x, y));
			}
			outputFile << '\n';
		}

		outputFile << "\nSolution time: " <<
			double(stop.QuadPart - start.QuadPart) / frequency.QuadPart / runsCount <<
			" seconds" << endl;

		break;
	} 

	return 0;
}
