#ifndef qtCWENO_nodeTag_H //include guard
#define qtCWENO_nodeTag_H

#include "main.h"
class nodeTag
{
    quadTreeId tree_ = null;
    int depth_ = 0;
    treeNodeId id_ = null;

public:
    nodeTag() {}; //default constructor
    nodeTag(quadTreeId t, int d, treeNodeId i) : tree_(t), depth_(d), id_(i) {}; //construct by all three values

    quadTreeId tree() const;
    int depth() const;
    treeNodeId id() const;
    bool operator == (const nodeTag& tag); //equality overload
    friend std::ostream& operator<<(std::ostream& os, const nodeTag& tag); //output overload
};

#endif
