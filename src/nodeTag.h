#ifndef qtCWENO_nodeTag_H //include guard
#define qtCWENO_nodeTag_H

#include "main.h"
class nodeTag
{
    treeId tree = null;
    int depth = 0;
    nodeId id = null;

public:
    nodeTag(); //default constructor
    nodeTag(treeId _tree_id, int _depth, nodeId _id); //construct by all three values
    treeId getTree() const;
    int getDepth() const;
    nodeId getId() const;
    bool operator == (const nodeTag& tag); //equality overload
    friend std::ostream& operator<<(std::ostream& os, const nodeTag& tag); //output overload
};

#endif
