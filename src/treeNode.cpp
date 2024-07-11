#include "main.h"
#include "cellBox.h"
#include "cellData.h"
#include "nodeEdge.h"
#include "nodeTag.h"
#include "quadTree.h"
#include "treeNode.h"

#include "Eigen\Dense"
#include <iomanip>

nodeTag treeNode::tag() const { return tag_; } //получение тэга 
void treeNode::setTag(nodeTag t) { tag_ = t; } //задание тэга
cellDataId treeNode::dataId() const { return dataId_; } //получение dataId
void treeNode::setDataId(cellDataId id) { dataId_ = id; } //задание dataId
cellBox treeNode::box() const {	return box_; } //получение box'a
void treeNode::setBox(cellBox b) { box_ = b; } //задание box'а
nodeTag treeNode::getNeighbour(Neighbour n) const {	return neighbours[static_cast<int>(n)]; }
void treeNode::setNeighbour(Neighbour n, nodeTag t) { neighbours[static_cast<int>(n)] = t; }


