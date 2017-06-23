// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/OPTICS-Clustering
// Distributed under the MIT Software License (X11 license).
// (See accompanying file LICENSE)


#pragma once

namespace optics {

template<typename T>
class Tree; //forward declaration of class Tree


template<typename T>
class Node {

	friend class Tree<T>;

public:
	Node( T data_ ) : data( data_ ) {}
	Node( T data_, std::vector<Node<T>> children_ ) : data( data_ ), children(children_) {}

	//void SetData( T data_ ) { data( data_ ); }
	T get_data() { return data; }
	T get_data() const { return data; }

	std::vector<Node<T>>& get_children() {
		return children;
	}
	const std::vector<Node<T>>& get_children() const{
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
	void append_subtree( std::vector<T>& result ) const{
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
	Tree() :root( { 0,0 } ) {}
	Tree( Node<T> root_ ) : root( root_ ) {}

	std::vector<T> flatten() {
		std::vector<T> result;
		root.append_subtree( result );
		return result;
	}

	Node<T>& get_root() {
		return root;
	}
	Node<T> get_root() const{
		return root;
	}

private:
	Node<T> root;

};


template<typename T>
std::size_t tree_depth( const Node<T>& root ) {
	std::size_t depth = 1;
	auto child_depths = fplus::transform( tree_depth<T>, root.get_children() );
	return 1 + (child_depths.empty() ? 0 : fplus::maximum( child_depths ));
			  
}


}//namespace optics