#include "main.h"
#include "NodeTag.h"


quadTreeId NodeTag::tree() const { return tree_; }
int NodeTag::depth() const { return depth_; }
treeNodeId NodeTag::id() const { return id_; }
bool NodeTag::isNull() const { return (tree() == null || id() == null); }

bool NodeTag::operator==(const NodeTag& rhs) //equality overload
{
    if (tree_ == rhs.tree_ && depth_ == rhs.depth_ && id_ == rhs.id_)
        return true;
    return false;
}

std::ostream& operator<<(std::ostream& os, const NodeTag& tag) //output overload
{
    os << "(" << tag.tree_ << ", " << tag.depth_ << ", " << tag.id_ << ")";
    return os;
}
