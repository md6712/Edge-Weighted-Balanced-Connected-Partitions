#include <iostream>
#include <algorithm>
using namespace std;

// Object class (data stored in the tree)
struct MyObject {
    int index;       // Unique index (key)
    void* item;      // Pointer to generic item

    MyObject(int idx, void* itm) : index(idx), item(itm) {}
};

// AVL Tree Node
class AVLTree {
private:
    struct Node {
        MyObject data;   // Object stored in the node
        Node* left;
        Node* right;
        int height;

        Node(MyObject obj) : data(obj), left(nullptr), right(nullptr), height(1) {}
    };

    Node* root;

    // Recursive function to delete the tree
    void deleteTree(Node* node) {
        if (!node) return;
        deleteTree(node->left);
        deleteTree(node->right);
        delete node; // Delete the current node
    }


    // Helper to get the height of a node
    int height(Node* n) {
        return n ? n->height : 0;
    }

    // Get the balance factor
    int balanceFactor(Node* n) {
        return n ? height(n->left) - height(n->right) : 0;
    }

    // Right rotate
    Node* rotateRight(Node* y) {
        Node* x = y->left;
        Node* T = x->right;
        x->right = y;
        y->left = T;
        y->height = max(height(y->left), height(y->right)) + 1;
        x->height = max(height(x->left), height(x->right)) + 1;
        return x;
    }

    // Left rotate
    Node* rotateLeft(Node* x) {
        Node* y = x->right;
        Node* T = y->left;
        y->left = x;
        x->right = T;
        x->height = max(height(x->left), height(x->right)) + 1;
        y->height = max(height(y->left), height(y->right)) + 1;
        return y;
    }

    // Insert operation
    Node* insert(Node* node, MyObject obj) {
        if (!node) return new Node(obj);

        if (obj.index < node->data.index)
            node->left = insert(node->left, obj);
        else if (obj.index > node->data.index)
            node->right = insert(node->right, obj);
        else
            return node; // Duplicate indexes are not allowed

        node->height = 1 + max(height(node->left), height(node->right));

        int balance = balanceFactor(node);

        // Balancing cases
        // Left-Left
        if (balance > 1 && obj.index < node->left->data.index)
            return rotateRight(node);

        // Right-Right
        if (balance < -1 && obj.index > node->right->data.index)
            return rotateLeft(node);

        // Left-Right
        if (balance > 1 && obj.index > node->left->data.index) {
            node->left = rotateLeft(node->left);
            return rotateRight(node);
        }

        // Right-Left
        if (balance < -1 && obj.index < node->right->data.index) {
            node->right = rotateRight(node->right);
            return rotateLeft(node);
        }

        return node;
    }

    // Inorder traversal
    void inorder(Node* node) {
        if (node) {
            inorder(node->left);
            cout << "Index: " << node->data.index << ", Item Address: " << node->data.item << endl;
            inorder(node->right);
        }
    }

    // Search for an object by index
    MyObject* search(Node* node, int index) {
        if (!node) return nullptr;
        if (index == node->data.index) return &(node->data);
        if (index < node->data.index) return search(node->left, index);
        return search(node->right, index);
    }

public:
    // Constructor
    AVLTree() : root(nullptr) {}

    // Destructor to clean up memory
    ~AVLTree() {
        deleteTree(root);
    }

    // Public insert method
    void insert(int index, void* item) {
        MyObject obj(index, item);
        root = insert(root, obj);
    }

    // Public inorder traversal method
    void inorder() {
        inorder(root);
    }

    // Public search method
    MyObject* search(int index) {
        return search(root, index);
    }
};


