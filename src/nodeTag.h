#ifndef qtCWENO_NodeTag_H //include guard
#define qtCWENO_NodeTag_H

#include "main.h"
class NodeTag
{
    quadTreeId tree_ = null;
    int depth_ = 0;
    treeNodeId id_ = null;

public:
    NodeTag() {}; //default constructor
    explicit NodeTag(quadTreeId t, int d, treeNodeId i) : tree_(t), depth_(d), id_(i) {}; //construct by all three values

    quadTreeId tree() const;
    int depth() const;
    treeNodeId id() const;
    bool isNull() const;
    bool operator==(const NodeTag& tag); //equality overload
    friend std::ostream& operator<<(std::ostream& os, const NodeTag& tag); //output overload
};

#endif
