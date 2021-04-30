#ifndef GUARD_QUEUE_H
#define GUARD_QUEUE_H


template<typename T>
class Queue {
public:
    Queue() : first(nullptr), last(nullptr), N(0) {}

    int size() const { return N; }
    bool isEmpty() const { return first == nullptr; }

    void enqueue(T val);
    T dequeue();

    // iterator
    class Iterator;
    Iterator begin();
    Iterator end();

private:
    // single Node in a list
    class Node {
    public:
        Node() : prev(nullptr), succ(nullptr) {};
        Node* prev;
        Node* succ;
        T val;
    };


    // first and last node in the list
    Node* first;
    Node* last;
    int N; // number of nodes
};


template<typename T>
void Queue<T>::enqueue(T val)
{
    Node* new_node = new Node();

    new_node->val = val;
    new_node->prev = last;
    if( last ) last->succ = new_node;
    last = new_node;
    if( !first ) first = last;
    ++N; 
}

template<typename T>
T Queue<T>::dequeue()
{
    T val = first->val;
    first = first->succ;
    --N;
    if(N == 0) last = nullptr;
    return val;

}

template<typename T>
class Queue<T>::Iterator {
private:
    Queue<T>::Node* current;

public:
    Iterator(Queue<T>::Node* p): current(p) {}

    Iterator& operator++() { current = current->succ; return *this;}
    T& operator*() const { return current->val; }
    bool hasNext() const { return current != nullptr; }

    bool operator==(const Iterator& b) const { return current == b.current; }
    bool operator!=(const Iterator& b) const { return current != b.current; }
};

template<class T>
typename Queue<T>::Iterator Queue<T>::begin()
{ return Iterator(first); }

template< class T>
typename Queue<T>::Iterator Queue<T>::end()
{ return Iterator(last); }

#endif
