#include "main.h"
#include "nodeTag.h"

nodeTag::nodeTag() {} //все дефолтные значени€ заданы в объ€влении структуры
nodeTag::nodeTag(treeId _tree, int _depth, nodeId _id) //construct by all three values 
{
    tree = _tree;
    depth = _depth;
    id = _id;
}
treeId nodeTag::getTree() const { return tree; }
int nodeTag::getDepth() const { return depth; }
nodeId nodeTag::getId() const { return id; }

bool nodeTag::operator==(const nodeTag& tag) //equality overload
{
    if (tree == tag.tree && depth == tag.depth && id == tag.id)
        return true;
    return false;
}

std::ostream& operator<<(std::ostream& os, const nodeTag& tag) //output overload
{
    os << "(" << tag.tree << ", " << tag.depth << ", " << tag.id << ")";
    return os;
}
