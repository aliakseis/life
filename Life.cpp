// Life.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <hash_set>
#include <iostream>

using stdext::hash_set;
using stdext::hash_compare;

using std::cout;


class Node;

struct IsNodePtrLess
{
	bool operator() (Node*, Node*) const;
};


namespace stdext 
{

size_t hash_value(Node*);

}

/**
*   This class contains the tree maintenance functions for quadtrees.
*/
class Node 
{

	typedef hash_set<Node*, hash_compare<Node*, IsNodePtrLess> > CacheType;

public:
	/**
	*   Construct a leaf cell.
	*/
	Node(bool living) 
	{
		nw = ne = sw = se = 0 ;
		level = 0 ;
		alive = living ;
		population = alive ? 1 : 0 ;
		result = 0;
	}
	/**
	*   Construct a node given four children.
	*/
	Node(Node* nw_, Node* ne_, Node* sw_, Node* se_) 
	{
		nw = nw_ ;
		ne = ne_ ;
		sw = sw_ ;
		se = se_ ;
		level = nw_->level + 1 ;
		population = nw->population + ne->population +
			sw->population + se->population ;
		alive = population > 0 ;
		result = 0;
	}
	/**
	*   Set a bit in this node in its relative coordinate system;
	*   returns a whole new node since our nodes are immutable.
	*
	*   In the recursive call, we simply adjust the coordinate system
	*   and call down a level.
	*/
	Node* setBit(int x, int y) {
		if (level == 0)
			//return new Node(true) ;
			return &aliveNode;
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
		if (level == 0)
			return alive ? 1 : 0 ;
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
	Node* emptyTree(int lev) {
		if (lev == 0)
			return create(false) ;
		Node* n = emptyTree(lev-1) ;
		return create(n, n, n, n) ;
	}
	/**
	*   Expand the universe; return a new node up one level with the
	*   current node in the center.  Requires us to disassemble the
	*   current node.
	*/
	Node* expandUniverse() {
		Node* border = emptyTree(level-1) ;
		return create(create(border, border, border, nw),
			create(border, border, ne, border),
			create(border, sw, border, border),
			create(se, border, border, border)) ;
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

	Node* sameGeneration() {
		return create(nw->se,
					  ne->sw,
					  sw->ne,
					  se->nw);
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
	Node* nextGeneration() {

		//if (level == 3)
		//	return create(nw->se,
		//				  ne->sw,
		//				  sw->ne,
		//				  se->nw);

		if (result != 0)
			return result ;
		if (population == 0)
			return result = nw ;
		if (level == 2)
			return result = slowSimulation() ;
		Node *n00 = nw->nextGeneration(),
			*n01 = horizontalForward(nw, ne),
			*n02 = ne->nextGeneration(),
			*n10 = verticalForward(nw, sw),
			*n11 = centerForward(),
			*n12 = verticalForward(ne, se),
			*n20 = sw->nextGeneration(),
			*n21 = horizontalForward(sw, se),
			*n22 = se->nextGeneration() ;

		//if (level == 3)
		//	return result = create(create(n00, n01, n10, n11)->sameGeneration(),
		//		create(n01, n02, n11, n12)->sameGeneration(),
		//		create(n10, n11, n20, n21)->sameGeneration(),
		//		create(n11, n12, n21, n22)->sameGeneration()) ;

		return result = create(create(n00, n01, n10, n11)->nextGeneration(),
			create(n01, n02, n11, n12)->nextGeneration(),
			create(n10, n11, n20, n21)->nextGeneration(),
			create(n11, n12, n21, n22)->nextGeneration()) ;
	}

	/**
	*   create functions.
	*/
	static Node* create(bool living) {
		return Node(living).intern() ;
	}
	static Node* create(Node* nw, Node* ne, Node* sw, Node* se) {
		return Node(nw, ne, sw, se).intern() ;
	}
	static Node* create() {
		return Node(false).emptyTree(3) ;
	}

	/**
	*   We need to provide a hashCode() and an equals() method to be
	*   able to hash these objects.
	*/

	size_t hashCode() 
	{
		if (level == 0)
			return population ;
		return size_t(nw) +
			11 * size_t(ne) +
			101 * size_t(sw) +
			1007 * size_t(se) ;
	}
	/*
	public bool equals(Object o) {
		TreeNode t = (TreeNode)o ;
		if (level != t.level)
			return false ;
		if (level == 0)
			return alive == t.alive ;
		return nw == t.nw && ne == t.ne && sw == t.sw && se == t.se ;
	}
	*/
	bool operator < (const Node& other) const
	{
		if (level < other.level)
			return true;
		if (level == 0)
			return alive < other.alive;

		return nw < other.nw 
			|| nw == other.nw && (ne < other.ne 
			|| ne == other.ne && (sw < other.sw || sw == other.sw && se < other.se));
	}

	/**
	*   Given a node, return the canonical one if it exists, or make it
	*   the canonical one.
	*/
	Node* intern() {

		CacheType::iterator it = cache.find(this);
		if (it != cache.end())
			return *it;

		Node* canon = new Node(*this);
		cache.insert(canon);

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
	Node* oneGen(int bitmask) {
		if (bitmask == 0)
			return create(false) ;
		int self = (bitmask >> 5) & 1 ;
		bitmask &= 0x757 ; // mask out bits we don't care about
		int neighborCount = 0 ;
		while (bitmask != 0) {
			neighborCount++ ;
			bitmask &= bitmask - 1 ; // clear least significant bit
		}
		if (neighborCount == 3 || (neighborCount == 2 && self != 0))
			return create(true) ;
		else
			return create(false) ;
	}
	/**
	*   At level 2, we can use slow simulation to compute the next
	*   generation.  We use bitmask tricks.
	*/
	Node* slowSimulation() {
		int allbits = 0 ;
		for (int y=-2; y<2; y++)
			for (int x=-2; x<2; x++)
				allbits = (allbits << 1) + getBit(x, y) ;
		return create(oneGen(allbits>>5), oneGen(allbits>>4),
			oneGen(allbits>>1), oneGen(allbits)) ;
	}

//private:
	Node *nw, *ne, *sw, *se ; // our children
	int level ;           // distance to root
	bool alive ;       // if leaf node, are we alive or dead?
	unsigned long population ;   // we cache the population here
	Node* result;

	static CacheType cache;
	static Node aliveNode;
};

Node::CacheType Node::cache;

Node Node::aliveNode(true);


bool IsNodePtrLess::operator() (Node* left, Node* right) const
{
	return *left < *right;
}


size_t stdext::hash_value(Node* pNode)
{
	return  pNode->hashCode();
}


class Universe
{
public:
	Universe()
	{
		generationCount = 0;
		root = Node::create();
	}

	void runStep() 
	{
		while (root->level < 3 ||
			root->nw->population != root->nw->se->se->population ||
			root->ne->population != root->ne->sw->sw->population ||
			root->sw->population != root->sw->ne->ne->population ||
			root->se->population != root->se->nw->nw->population)
			root = root->expandUniverse() ;

		unsigned long stepSize = 1ul << (root->level - 2);

		root = root->nextGeneration() ;
		generationCount += stepSize ;
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

//private:
	Node* root;
	unsigned long generationCount;
};


int main(int argc, char* argv[])
{
	Universe universe;

	universe.setBit(1, 0);
	universe.setBit(2, 0);
	universe.setBit(0, 1);
	universe.setBit(1, 1);
	universe.setBit(1, 2);

	universe.root = universe.root->expandUniverse() ;
	universe.root = universe.root->expandUniverse() ;
	universe.root = universe.root->expandUniverse() ;
	universe.root = universe.root->expandUniverse() ;

	//do
	{
		universe.runStep();
		cout << universe.generationCount << ' ' << universe.root->population << '\n';
	}
	//while (universe.generationCount < 1000);

//*
	int radius = 1 << (universe.root->level - 1);
	for (int y = -radius; y < radius; ++y)
	{
		for (int x = -radius; x < radius; ++x)
		{
			cout << (universe.root->getBit(x, y)? '*' : '.');
		}
		cout << '\n';
	}
//*/
	return 0;
}

