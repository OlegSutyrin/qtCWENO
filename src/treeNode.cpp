#include "main.h"
#include "cellBox.h"
#include "cellData.h"
#include "nodeEdge.h"
#include "nodeTag.h"
#include "quadTree.h"
#include "treeNode.h"

#include "Eigen\Dense"
#include <iomanip>

nodeTag treeNode::tag() const { return tag_; } //��������� ���� 
void treeNode::setTag(nodeTag t) { tag_ = t; } //������� ����
cellDataId treeNode::dataId() const { return dataId_; } //��������� dataId
void treeNode::setDataId(cellDataId id) { dataId_ = id; } //������� dataId
cellBox treeNode::box() const {	return box_; } //��������� box'a
void treeNode::setBox(cellBox b) { box_ = b; } //������� box'�
nodeTag treeNode::getNeighbour(Neighbour n) const {	return neighbours[static_cast<int>(n)]; }
void treeNode::setNeighbour(Neighbour n, nodeTag t) { neighbours[static_cast<int>(n)] = t; }


