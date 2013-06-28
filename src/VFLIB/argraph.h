/*------------------------------------------------------------------
 * argraph.h
 * Interface of argraph.cc
 * Definition of a class representing an Attributed Relational 
 * Graph (ARG).
 *
 * Author: P. Foggia
 *-----------------------------------------------------------------*/



/*-----------------------------------------------------------------
 *   IMPLEMENTATION NOTES
 *
 * There is an implementation class, ARGraph_impl, which does
 * the real work and represents the node and edge attributes
 * using void* pointers; 
 * A template-based interface class ARGraph<NodeType,EdgeType>
 * provides type checking. Attributes are always managed 
 * via pointers. 
 *
 * The user has two choices for the attribute allocation strategy:
 * either the attributes are not owned by the graph, and someone
 * else has to provide for their deallocation, or the attributes
 * are owned by the graph, and are deallocated by means of a
 * user provided object or function.
 *
 * The graphs, once created, are immutable, in the sense that
 * graph-editing operations are not provided.
 * This allows the use of an internal representation particularly
 * suited for efficient graph matching, which is the primary target
 * of this program.
 *
 * An abstract class, ARGLoader, is defined to allow the
 * use of different file formats for loading the graphs.
 * The loader is queried by the ARGraph to acquire the information
 * needed for building the graph.
 * A simple ARGLoader based on iostream.h is provided in
 * argloader.h; see also argedit.h/argedit.cc for a class which can
 * be used as a base for a new ARGLoader.
 --------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 *   MORE IMPLEMENTATION DETAILS
 * Edges and edge attributes (pointers) are stored into sorted 
 * vectors associated to each node, and are looked for
 * using binary search. 
 * 
 * Nodes are identified using the type node_id, which is currently
 * unsigned short; the special value NULL_NODE is used as null
 * value for this type.
 *
 * Bound checks are performed using the assert macro. They can be
 * disabled by ensuring the macro NDEBUG is defined during 
 * compilation.
 *
 * NOTE: Differently from the previous versions  of this library 
 * (before version 2.0), there is no more an adjacency matrix to 
 * check for the existence of a node.
 --------------------------------------------------------------------*/

#ifndef ARGRAPH_H
#define ARGRAPH_H

#include <assert.h>
#include <stddef.h>

/*--------------------------------------------------------
 * General definitions
 -------------------------------------------------------*/
#ifdef NEED_BOOL
#ifndef BOOL_DEFINED
#define BOOL_DEFINED
typedef int bool;
const false=0;
const true=1
#endif
#endif

#ifndef BYTE_DEFINED
#define BYTE_DEFINED
typedef unsigned char byte;
#endif

typedef unsigned short node_id;
const node_id NULL_NODE=0xFFFF;

/*----------------------------------------------------------------
 * Abstract class ARGLoader. Allows to construct an ARGraph
 ---------------------------------------------------------------*/
class ARGLoader
  { public:
      virtual ~ARGLoader() {}


      virtual int NodeCount() = 0;
      virtual void *GetNodeAttr(node_id node) = 0;
      virtual int OutEdgeCount(node_id node) = 0;
      virtual node_id GetOutEdge(node_id node, int i, void **pattr)=0;
  };


/*------------------------------------------------------------
 * Abstract class AttrDestroyer. Allows to destroy
 * a node or edge attribute.
 -----------------------------------------------------------*/
class AttrDestroyer
  { public:
      virtual ~AttrDestroyer() {
    	  // Do nothing
      };
      virtual void destroy(void *attr)=0;
      // linh
      virtual AttrDestroyer* clone(){
    	  return NULL;
      }
      // end linh
  };

/*----------------------------------------------------------
 * Abstract class AttrComparator. Allows to check for
 * compatibility between two node or egde attributes.
 ---------------------------------------------------------*/
class AttrComparator
  { public:
      virtual ~AttrComparator() {
    	  // Do nothing
      };
      virtual bool compatible(void *attr1, void *attr2)=0;
      // linh
      virtual AttrComparator* clone() {
    	  return NULL;
      }
      // end linh
  };

/*--------------------------------------------------------
 * Declaration of class ARGraph_impl
 * This is the real representation of an ARG, but it is
 * intended to be used thru the interface ARGraph
 -------------------------------------------------------*/
class ARGraph_impl
  { public:
      typedef void *param_type;
      typedef void (*edge_visitor)(ARGraph_impl *g, 
                                   node_id n1, node_id n2, void *attr,
                                   param_type param);
      typedef void (*node_destroy_fn)(void *);
      typedef void (*edge_destroy_fn)(void *);

      typedef bool (*node_compat_fn)(void *, void *);
      typedef bool (*edge_compat_fn)(void *, void *);



    protected:
      typedef short count_type;

      int n;              /* number of nodes  */
      void* *attr;        /* node attributes  */
      count_type  *in_count;  /* number of 'in' edges for each node */
      node_id **in;       /* nodes connected by 'in' edges to each node */
      void* **in_attr;    /* Edge attributes for 'in' edges */
      count_type  *out_count; /* number of 'out edges for each node */
      node_id **out;      /* nodes connected by 'out' edges to each node */
      void* **out_attr;   /* Edge attributes for 'out' edges */

      AttrDestroyer *node_destroyer;  // Used to clean up node attrs
      AttrDestroyer *edge_destroyer;  // Used to clean up edge attrs
      AttrComparator *node_comparator; // Used to test node attr. compat.
      AttrComparator *edge_comparator; // Used to test edge attr. compat.

    public: // was protected:
      void SetNodeDestroyer(AttrDestroyer *);
      void SetEdgeDestroyer(AttrDestroyer *);
      void SetNodeDestroy(node_destroy_fn fn); // For older versions
      void SetEdgeDestroy(edge_destroy_fn fn); // For older versions

      void SetNodeComparator(AttrComparator *);
      void SetEdgeComparator(AttrComparator *);
      void SetNodeCompat(node_compat_fn); // For older versions
      void SetEdgeCompat(edge_compat_fn); // For older versions

    public:
      // linh
      ARGraph_impl();
      void setLoader(ARGLoader *loader);
      void ClearAll();
      // end linh

      ARGraph_impl(ARGLoader *loader);
      virtual ~ARGraph_impl();


      int NodeCount();

      void *GetNodeAttr(node_id i);
      void SetNodeAttr(node_id i, void *attr, bool destroyOld=false);

      bool HasEdge(node_id n1, node_id n2);
      bool HasEdge(node_id n1, node_id n2, void **pattr);
      void *GetEdgeAttr(node_id n1, node_id n2);
      void SetEdgeAttr(node_id n1, node_id n2, void *attr, 
                       bool destroyOld=false);

      int InEdgeCount(node_id node);
      int OutEdgeCount(node_id node); 
      int EdgeCount(node_id node); 
      node_id GetInEdge(node_id node, int i);
      node_id GetInEdge(node_id node, int i, void **pattr);
      node_id GetOutEdge(node_id node, int i);
      node_id GetOutEdge(node_id node, int i, void **pattr);


      void VisitInEdges(node_id node, edge_visitor vis, param_type param);
      void VisitOutEdges(node_id node, edge_visitor vis, param_type param);
      void VisitEdges(node_id node, edge_visitor vis, param_type param);

      virtual bool CompatibleNode(void *attr1, void *attr2);
      virtual bool CompatibleEdge(void *attr1, void *attr2);

    protected:
      virtual void DestroyNode(void *attr);
      virtual void DestroyEdge(void *attr);
  };

/*
 * The following typedef is for compatibility
 */
typedef ARGraph_impl Graph;


/*--------------------------------------------------------
 * Declaration of class ARGraph
 * This is the interface to manage ARGs.
 -------------------------------------------------------*/
template <class Node, class Edge>
class ARGraph: public ARGraph_impl
  { public:
      typedef Node node_type;
      typedef Edge edge_type;

      typedef void *param_type;
      typedef void (*edge_visitor)(ARGraph *g, 
                                   node_id n1, node_id n2, Edge *attr,
                                   param_type param);
      typedef void (*node_destroy_fn)(Node *);
      typedef void (*edge_destroy_fn)(Edge *);

      typedef bool (*node_compat_fn)(Node *, Node *);
      typedef bool (*edge_compat_fn)(Edge *, Edge *);

      typedef ARGraph_impl impl;


      /*-------- methods ----------*/
      // linh
      ARGraph(): impl() {};
      ARGraph(const ARGraph& another_argraph);
      ARGraph& operator = (const ARGraph& another_argraph);
      // end linh

      ARGraph(ARGLoader *loader) : impl(loader) {}

      void SetNodeDestroy(node_destroy_fn fn)
          { ARGraph_impl::SetNodeDestroy((ARGraph_impl::node_destroy_fn) fn); }
      void SetEdgeDestroy(edge_destroy_fn fn)
          { ARGraph_impl::SetEdgeDestroy((ARGraph_impl::edge_destroy_fn) fn); }

      void SetNodeCompat(node_compat_fn fn)
          { ARGraph_impl::SetNodeCompat((ARGraph_impl::node_compat_fn) fn); }
      void SetEdgeCompat(edge_compat_fn fn)
          { ARGraph_impl::SetEdgeCompat((ARGraph_impl::edge_compat_fn) fn); }


      bool HasEdge(node_id n1, node_id n2, Edge **pattr)
          { return ARGraph_impl::HasEdge(n1, n2, (void**)pattr); }

      Node *GetNodeAttr(node_id i)
          { return (Node*)ARGraph_impl::GetNodeAttr(i); }

      Edge *GetEdgeAttr(node_id n1, node_id n2)
                 { return (Edge*)ARGraph_impl::GetEdgeAttr(n1,n2); }

      // Linh
      node_id GetInEdge(node_id node, int i)
                 { return ARGraph_impl::GetInEdge(node, i); }
      node_id GetOutEdge(node_id node, int i)
                 { return ARGraph_impl::GetOutEdge(node, i); }
      // End Linh
      node_id GetInEdge(node_id node, int i, Edge **pattr)
                 { return ARGraph_impl::GetInEdge(node, i, (void**)pattr); }
      node_id GetOutEdge(node_id node, int i, Edge **pattr)
                 { return ARGraph_impl::GetOutEdge(node, i, (void**)pattr); }


      void VisitInEdges(node_id node, edge_visitor vis, param_type param)
         { ARGraph_impl::VisitInEdges(node, (ARGraph_impl::edge_visitor)vis, 
                                      param); }
      void VisitOutEdges(node_id node, edge_visitor vis, param_type param)
         { ARGraph_impl::VisitOutEdges(node, (ARGraph_impl::edge_visitor)vis, 
                                       param); }
      void VisitEdges(node_id node, edge_visitor vis, param_type param)
         { ARGraph_impl::VisitEdges(node, (ARGraph_impl::edge_visitor)vis, 
                                    param); }


      bool CompatibleNode(Node *attr1, Node *attr2)
         { return ARGraph_impl::CompatibleNode(attr1, attr2); }
      bool CompatibleEdge(Edge *attr1, Edge *attr2)
         { return ARGraph_impl::CompatibleEdge(attr1, attr2); }
  };

// linh
template <class Node, class Edge>
ARGraph<Node, Edge>::ARGraph(const ARGraph<Node, Edge>& another_argraph) {

	n = another_argraph.n;
	attr = (void**) new Node*[n];
	for (int i = 0; i < n; i++)
		// Linh
		attr[i] = (void*) new Node(*((Node*)another_argraph.attr[i]));

	out_count = new count_type[n];
	in_count = new count_type[n];
	out = new node_id*[n];
	in = new node_id*[n];
	out_attr = new void**[n];
	in_attr = new void**[n];
	for (int i = 0; i < n; i++){
		out_count[i] = another_argraph.out_count[i];
		out[i] = new node_id[out_count[i]];
		out_attr[i] = new void*[out_count[i]];
		for (int j = 0; j < out_count[i]; j++) {
			out[i][j] = another_argraph.out[i][j];
			out_attr[i][j] = (void*) new Edge(*((Edge*)another_argraph.out_attr[i][j]));
		}
		in_count[i] = another_argraph.in_count[i];
		in[i] = new node_id[in_count[i]];
		in_attr[i] = new void*[in_count[i]];
		for (int j = 0; j < in_count[i]; j++) {
			in[i][j] = another_argraph.in[i][j];
			in_attr[i][j] = (void*) new Edge(*((Edge*)another_argraph.in_attr[i][j]));
		}
	}
	node_destroyer = another_argraph.node_destroyer->clone();
	edge_destroyer = another_argraph.edge_destroyer->clone();
	node_comparator = another_argraph.node_comparator->clone();
	edge_comparator = another_argraph.edge_comparator->clone();
}

template <class Node, class Edge>
ARGraph<Node, Edge>& ARGraph<Node, Edge>::operator = (const ARGraph<Node, Edge>& another_argraph) {

	if (this != &another_argraph) {
		ClearAll();	// Deallocate memory
		n = another_argraph.n;
		attr = (void**) new Node*[n];
		for (int i = 0; i < n; i++)
			// Linh
			attr[i] = (void*) new Node(*((Node*)another_argraph.attr[i]));

		out_count = new count_type[n];
		in_count = new count_type[n];
		out = new node_id*[n];
		in = new node_id*[n];
		out_attr = new void**[n];
		in_attr = new void**[n];
		for (int i = 0; i < n; i++){
			out_count[i] = another_argraph.out_count[i];
			out[i] = new node_id[out_count[i]];
			out_attr[i] = new void*[out_count[i]];
			for (int j = 0; j < out_count[i]; j++) {
				out[i][j] = another_argraph.out[i][j];
				out_attr[i][j] = (void*) new Edge(*((Edge*)another_argraph.out_attr[i][j]));
			}
			in_count[i] = another_argraph.in_count[i];
			in[i] = new node_id[in_count[i]];
			in_attr[i] = new void*[in_count[i]];
			for (int j = 0; j < in_count[i]; j++) {
				in[i][j] = another_argraph.in[i][j];
				in_attr[i][j] = (void*) new Edge(*((Edge*)another_argraph.in_attr[i][j]));
			}
		}
		node_destroyer = another_argraph.node_destroyer->clone();
		edge_destroyer = another_argraph.edge_destroyer->clone();
		node_comparator = another_argraph.node_comparator->clone();
		edge_comparator = another_argraph.edge_comparator->clone();
	}
	return *this;
}
// end linh

/*---------------------------------------------------------------------
 * Declaration of classes FunctionAttrDestroyer and
 * FunctionAttrComparator, used to create an AttrDestroyer or
 * an AttrComparator from a C function.
 --------------------------------------------------------------------*/
class FunctionAttrDestroyer: public AttrDestroyer
  { private:
      void (*func)(void *);
    public:
      FunctionAttrDestroyer(void (*fn)(void *));
      virtual void destroy(void *attr);
  };

class FunctionAttrComparator: public AttrComparator
  { private:
      bool (*func)(void *, void *);
    public:
      FunctionAttrComparator(bool (*fn)(void *, void *));
      virtual bool compatible(void *attr1, void *attr2);
  };


/*---------------------------------------------------------------------
 * INLINE METHODS
 ---------------------------------------------------------------------*/

/*-----------------------------------------------
 * Class ARGraph_impl
 ----------------------------------------------*/

/*-----------------------------------------------
 * Number of nodes in the graph
 ----------------------------------------------*/
inline int ARGraph_impl::NodeCount() 
  { return n;
  } 

/*----------------------------------------------
 * Attribute of a node
 ---------------------------------------------*/
inline void * ARGraph_impl::GetNodeAttr(node_id i) 
  { assert(i<n);
    return attr[i];
  }

/*----------------------------------------------
 * Check the presence of an edge
 ---------------------------------------------*/
inline bool ARGraph_impl::HasEdge(node_id n1, node_id n2)
    { return HasEdge(n1, n2, NULL);
    }

/*----------------------------------------------
 * Gets the attribute of an edge
 ---------------------------------------------*/
inline void* ARGraph_impl::GetEdgeAttr(node_id n1, node_id n2)
    { void *attr;
      if (HasEdge(n1, n2, &attr))
         return attr;
      else
         return NULL;
    }


/*------------------------------------------------
 * Number of edges going into a node
 ------------------------------------------------*/
inline int ARGraph_impl::InEdgeCount(node_id node) 
  { assert(node<n); 
    return in_count[node]; 
  }


/*------------------------------------------------
 * Number of edges departing from a node
 ------------------------------------------------*/
inline int ARGraph_impl::OutEdgeCount(node_id node) 
  { assert(node<n); 
    return out_count[node]; 
  }


/*-------------------------------------------------
 * Number of edges touching a node
 ------------------------------------------------*/
inline int ARGraph_impl::EdgeCount(node_id node) 
  { assert(node<n); 
    return in_count[node]+out_count[node]; 
  }

/*------------------------------------------------
 * Gets the other end of an edge entering a node
 -----------------------------------------------*/
inline node_id ARGraph_impl::GetInEdge(node_id node, int i)
  { assert(node<n);
    assert(i<in_count[node]);
    return in[node][i];
  }

/*------------------------------------------------
 * Gets the other end of an edge entering a node
 * Also gets the attribute of the edge
 -----------------------------------------------*/
inline node_id ARGraph_impl::GetInEdge(node_id node, int i, 
                                       void **pattr)
  { assert(node<n);
    assert(i<in_count[node]);
    *pattr = in_attr[node][i];
    return in[node][i];
  }


/*------------------------------------------------
 * Gets the other end of an edge leaving a node
 -----------------------------------------------*/
inline node_id ARGraph_impl::GetOutEdge(node_id node, int i)
  { assert(node<n);
    assert(i<out_count[node]);
    return out[node][i];
  }

/*------------------------------------------------
 * Gets the other end of an edge leaving a node
 * Also gets the attribute of the edge
 -----------------------------------------------*/
inline node_id ARGraph_impl::GetOutEdge(node_id node, int i, 
                                        void **pattr)
  { assert(node<n);
    assert(i<out_count[node]);
    *pattr = out_attr[node][i];
    return out[node][i];
  }

/*-----------------------------------------------------------
 * Checks if two node attributes are compatible
 ----------------------------------------------------------*/
inline bool ARGraph_impl::CompatibleNode(void *attr1, void *attr2)
  { if (node_comparator==NULL)
      return true;
    else
      return node_comparator->compatible(attr1, attr2);
  }


/*-----------------------------------------------------------
 * Checks if two edge attributes are compatible
 ----------------------------------------------------------*/
inline bool ARGraph_impl::CompatibleEdge(void *attr1, void *attr2)
  { if (edge_comparator==NULL)
      return true;
    else
      return edge_comparator->compatible(attr1, attr2);
  }


#endif
/* defined ARGRAPH_H */


