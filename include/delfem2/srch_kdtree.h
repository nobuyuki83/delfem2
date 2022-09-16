#ifndef DFM2_KDT_H
#define DFM2_KDT_H

#include <cassert>
#include <vector>
#include <algorithm>

namespace delfem2 {

template <typename VEC, typename SCALAR = typename VEC::Scalar>
struct KDTNodePoint {
	KDTNodePoint* parent	 = nullptr;
	KDTNodePoint* rightChild = nullptr;
	KDTNodePoint* leftChild	 = nullptr;
	VEC separator;
};

template <typename VEC, typename SCALAR = typename VEC::Scalar>
void ConstructNode(KDTNodePoint<VEC>* currentRoot, std::vector<VEC> pointList, const int beginInd, const int endInd, const int depth = 0)
{

	if ((endInd - beginInd) == 1) {
		currentRoot->separator = pointList[beginInd];
		return;
	}

	//sort
	{
		auto begin = pointList.begin();
		std::advance(begin, beginInd);

		auto end = pointList.begin();
		std::advance(end, endInd);

		if (depth % 3 == 0) {
			std::sort(begin, end,
			    [](VEC a, VEC b) -> bool { return a[0] > b[0]; });
		} else if (depth % 3 == 1) {
			std::sort(begin, end,
			    [](VEC a, VEC b) -> bool { return a[1] > b[1]; });
		} else {
			std::sort(begin, end,
			    [](VEC a, VEC b) -> bool { return a[2] > b[2]; });
		}
	}

	int midInd	       = (endInd - beginInd) / 2 + beginInd;
	currentRoot->separator = pointList[midInd];

	if (beginInd != midInd) {
		currentRoot->rightChild		= new KDTNodePoint<VEC>;
		currentRoot->rightChild->parent = currentRoot;
		ConstructNode(currentRoot->rightChild, pointList, beginInd, midInd, depth + 1);
	}

	if ((midInd + 1) != endInd) {
		currentRoot->leftChild	       = new KDTNodePoint<VEC>;
		currentRoot->leftChild->parent = currentRoot;
		ConstructNode(currentRoot->leftChild, pointList, midInd + 1, endInd, depth + 1);
	}
}

template <typename VEC, typename SCALAR = typename VEC::Scalar>
KDTNodePoint<VEC>* ConstructKDT(const std::vector<VEC>& Points)
{
	assert(!Points.empty());

	KDTNodePoint<VEC>* root = new KDTNodePoint<VEC>;
	root->parent		= nullptr;

	std::vector<VEC> pointList = Points;

	ConstructNode<VEC>(root, pointList, 0, pointList.size());

	return root;
}

template <typename VEC, typename SCALAR = typename VEC::Scalar>
VEC SearchNearstPoint(const KDTNodePoint<VEC>* const root, const VEC Q, const int depth, double& closestdist)
{
	double localclosestdist = closestdist;
	VEC retv;

	const KDTNodePoint<VEC>* nextNode    = nullptr;
	const KDTNodePoint<VEC>* siblingNode = nullptr;
	double distP;
	if (depth % 3 == 0)
		distP = Q[0] - root->separator[0];
	else if (depth % 3 == 1)
		distP = Q[1] - root->separator[1];
	else
		distP = Q[2] - root->separator[2];

	if (root->rightChild != nullptr && root->leftChild != nullptr) {
		if (distP > 0.0) {
			nextNode    = root->rightChild;
			siblingNode = root->leftChild;
		} else {
			nextNode    = root->leftChild;
			siblingNode = root->rightChild;
		}
	} else if (root->rightChild != nullptr) {
		nextNode    = root->rightChild;
		siblingNode = root->leftChild; //nullptr
	} else if (root->leftChild != nullptr) {
		nextNode    = root->leftChild;
		siblingNode = root->rightChild; //nullptr
	} else {
	}

	const double dist0 = (root->separator - Q).norm();
	if (dist0 < localclosestdist) {
		retv		 = root->separator;
		localclosestdist = dist0;
	}

	if (nextNode != nullptr) {
		double dist1	    = localclosestdist;
		const VEC closestP1 = SearchNearstPoint(nextNode, Q, depth + 1, dist1);
		if (dist1 > -0.1 && dist1 < localclosestdist) {
			retv		 = closestP1;
			localclosestdist = dist1;
		}
	}

	if (std::abs(distP) < localclosestdist && siblingNode != nullptr) {
		double dist2	    = localclosestdist;
		const VEC closestP2 = SearchNearstPoint(siblingNode, Q, depth + 1, dist2);
		if (dist2 > -0.1 && dist2 < localclosestdist) {
			retv		 = closestP2;
			localclosestdist = dist2;
		}
	}

	if (localclosestdist > closestdist) {
		closestdist = -1.0;
		return { 0.0, 0.0, 0.0 };
	} else {
		closestdist = localclosestdist;
		return retv;
	}
}

template <typename VEC, typename SCALAR = typename VEC::Scalar>
VEC SearchNearstPoint(const KDTNodePoint<VEC>* const root, const VEC Q)
{
	double temp = std::numeric_limits<double>::max();
	return SearchNearstPoint(root, Q, 0, temp);
}

}

#endif //DFM2_KDT_H
