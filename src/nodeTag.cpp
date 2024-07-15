#include "main.h"
#include "NodeTag.h"


quadTreeId NodeTag::tree() const { return tree_; }
int NodeTag::depth() const { return depth_; }
treeNodeId NodeTag::id() const { return id_; }
bool NodeTag::isNull() const { return (tree() == null || id() == null); }

bool NodeTag::operator==(const NodeTag& tag) //equality overload
{
    if (tree_ == tag.tree_ && depth_ == tag.depth_ && id_ == tag.id_)
        return true;
    return false;
}

std::ostream& operator<<(std::ostream& os, const NodeTag& tag) //output overload
{
    os << "(" << tag.tree_ << ", " << tag.depth_ << ", " << tag.id_ << ")";
    return os;
}
