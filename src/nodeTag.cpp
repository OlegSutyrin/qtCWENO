#include "main.h"
#include "nodeTag.h"


quadTreeId nodeTag::tree() const { return tree_; }
int nodeTag::depth() const { return depth_; }
treeNodeId nodeTag::id() const { return id_; }

bool nodeTag::operator==(const nodeTag& tag) //equality overload
{
    if (tree_ == tag.tree_ && depth_ == tag.depth_ && id_ == tag.id_)
        return true;
    return false;
}

std::ostream& operator<<(std::ostream& os, const nodeTag& tag) //output overload
{
    os << "(" << tag.tree_ << ", " << tag.depth_ << ", " << tag.id_ << ")";
    return os;
}
