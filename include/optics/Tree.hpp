// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/OPTICS-Clustering
// Distributed under the MIT Software License (X11 license).
// (See accompanying file LICENSE)


#pragma once

namespace optics {


template<typename T>
class Node {
	template<typename T>
	friend class Tree;
public:
	Node( T data_ ) : data( data_ ) {}

	//void SetData( T data_ ) { data( data_ ); }
	T get_data( ) { return data; }
	
	std::vector<Node<T>>& get_children() {
		return children;
	}

	void add_child( const Node<T>& new_child ) {
		children.push_back( new_child );
	}
	void add_children( const std::vector<Node<T>>& new_children ) {
		for ( const auto& c : new_children ) {
			children.push_back( c );
		}
	}

private:
	void Node::append_subtree( std::vector<T>& result ) const{
		result.push_back( data );
		for ( const auto& c : children ) {
			c.append_subtree(result);
		}
	};

	std::vector<Node<T>> children;
	T data;
};


template<typename T>
class Tree {
public:
	Tree( Node<T> root_ ) : root( root_ ) {}

	std::vector<T> Tree::flatten() {
		std::vector<T> result;
		root.append_subtree( result );
		return result;
	}

	Node<T>& get_root() {
		return root;
	}

private:
	Node<T> root;

};

}//namespace optics