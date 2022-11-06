#ifndef DFM2_KDT_TRI_H
#define DFM2_KDT_TRI_H

#include <cassert>
#include <vector>
#include <algorithm>
#include <cmath>

#include <iostream>

namespace delfem2 {

template <typename VEC, typename SCALAR = typename VEC::Scalar>
struct KDTTriangle {
    private:
	int root;

	struct KDTNodeTriangle {
		int TriInd;
		SCALAR AABBmax[3];
		SCALAR AABBmin[3];
		VEC tricenter;

		int parent	     = -1;
		int rightchild	     = -1;
		float realboundright = 0.0;
		int leftchild	     = -1;
		float realboundleft  = 0.0;

		KDTNodeTriangle(const int TriInd, const std::vector<VEC>& PointList, const std::vector<unsigned int>& TriList)
		    : TriInd(TriInd)
		    , AABBmax { std::max({ PointList[TriList[3 * TriInd + 0]][0], PointList[TriList[3 * TriInd + 1]][0], PointList[TriList[3 * TriInd + 2]][0] }),
			    std::max({ PointList[TriList[3 * TriInd + 0]][1], PointList[TriList[3 * TriInd + 1]][1], PointList[TriList[3 * TriInd + 2]][1] }),
			    std::max({ PointList[TriList[3 * TriInd + 0]][2], PointList[TriList[3 * TriInd + 1]][2], PointList[TriList[3 * TriInd + 2]][2] }) }
		    , AABBmin { std::min({ PointList[TriList[3 * TriInd + 0]][0], PointList[TriList[3 * TriInd + 1]][0], PointList[TriList[3 * TriInd + 2]][0] }),
			    std::min({ PointList[TriList[3 * TriInd + 0]][1], PointList[TriList[3 * TriInd + 1]][1], PointList[TriList[3 * TriInd + 2]][1] }),
			    std::min({ PointList[TriList[3 * TriInd + 0]][2], PointList[TriList[3 * TriInd + 1]][2], PointList[TriList[3 * TriInd + 2]][2] }) }
		    , tricenter((1.0 / 3.0) * (PointList[TriList[3 * TriInd + 0]] + PointList[TriList[3 * TriInd + 1]] + PointList[TriList[3 * TriInd + 2]]))
		{
		}

		auto distance(VEC Q)
		{
			struct {
				SCALAR max, min;
			} dist;

			dist.min = (Q - VEC { AABBmax[0], AABBmax[1], AABBmax[2] }).norm();
			dist.max = (Q - VEC { AABBmax[0], AABBmax[1], AABBmax[2] }).norm();

			for (int i = 0; i < 8; i++) {
				SCALAR hogex = ((i >> 0) & 1) ? AABBmax[0] : AABBmin[0];
				SCALAR hogey = ((i >> 1) & 1) ? AABBmax[1] : AABBmin[1];
				SCALAR hogez = ((i >> 2) & 1) ? AABBmax[2] : AABBmin[2];

				SCALAR hogedist = (Q - VEC { hogex, hogey, hogez }).norm();

				if (hogedist < dist.min)
					dist.min = hogedist;
				if (hogedist > dist.max)
					dist.max = hogedist;
			}

			if (Q[0] < AABBmax[0] && Q[0] > AABBmin[0] && Q[1] < AABBmax[1] && Q[1] > AABBmin[1] && Q[2] < AABBmax[2] && Q[2] > AABBmin[2])
				dist.min = 0.0;

			return dist;
		}
	};

	std::vector<KDTNodeTriangle> NodeList;

	struct iss {
		int index;
		SCALAR max;
		SCALAR min;
	};

	int SearchNearestRecursive(const int currentInd, const VEC Query, std::vector<iss>& Candidates, SCALAR& distmax, SCALAR& distmin, const unsigned int depth)
	{
		SCALAR localmax = distmax;
		SCALAR localmin = distmin;

		int onnodeInd  = -1;
		int outnodeInd = -1;

		SCALAR distfromPlane;
		if (depth % 3 == 0)
			distfromPlane = Query[0] - NodeList[currentInd].tricenter[0];
		else if (depth % 3 == 1)
			distfromPlane = Query[1] - NodeList[currentInd].tricenter[1];
		else
			distfromPlane = Query[2] - NodeList[currentInd].tricenter[2];

		SCALAR outrealbound;

		if (distfromPlane > 0.0) {
			onnodeInd    = NodeList[currentInd].rightchild;
			outnodeInd   = NodeList[currentInd].leftchild;
			outrealbound = NodeList[currentInd].realboundleft;
		} else {
			onnodeInd    = NodeList[currentInd].leftchild;
			outnodeInd   = NodeList[currentInd].rightchild;
			outrealbound = NodeList[currentInd].realboundright;
		}

		auto distfromcenter = NodeList[currentInd].distance(Query);
		if (distfromcenter.max < localmin) {
			localmax = distfromcenter.max;
			localmin = distfromcenter.min;
			Candidates.clear();
			Candidates.push_back({ currentInd, distfromcenter.max, distfromcenter.min });
		} else if (distfromcenter.min < localmin) {
			localmax = distfromcenter.max;
			localmin = distfromcenter.min;
			std::remove_if(Candidates.begin(), Candidates.end(), [localmax](iss x) {
				return localmax < x.min - 0.001;
			});
			Candidates.push_back({ currentInd, distfromcenter.max, distfromcenter.min });
		} else if (distfromcenter.min < localmax) {
			Candidates.push_back({ currentInd, distfromcenter.max, distfromcenter.min });
		}

		if (onnodeInd != -1) {
			SearchNearestRecursive(onnodeInd, Query, Candidates, localmax, localmin, depth + 1);
		}

		if (outnodeInd != -1) {
			SCALAR realdist;
			if (depth % 3 == 0)
				realdist = Query[0] - outrealbound;
			else if (depth % 3 == 1)
				realdist = Query[1] - outrealbound;
			else
				realdist = Query[2] - outrealbound;

			if ((realdist * distfromPlane) < 0.001 || std::abs(realdist) < localmax + 0.001) {
				SearchNearestRecursive(outnodeInd, Query, Candidates, localmax, localmin, depth + 1);
			}
		}

		distmax = localmax;
		distmin = localmin;

		return 1;
	}

    public:
	KDTTriangle(const std::vector<VEC>& PL, const std::vector<unsigned int>& TL)
	{

		assert(!PL.empty());
		assert(!TL.empty());

		for (int i = 0; i < TL.size() / 3; i++) {
			NodeList.emplace_back(i, PL, TL);
		}

		root = ConstructKDTree(0, NodeList.size(), 0);
	}

	int ConstructKDTree(const int beginInd, const int endInd, const int depth)
	{
		if ((endInd - beginInd) == 1) {
			return beginInd;
		}

		auto begin = NodeList.begin();
		std::advance(begin, beginInd);

		auto end = NodeList.begin();
		std::advance(end, endInd);

		if (depth % 3 == 0) {
			std::sort(begin, end, [](const KDTNodeTriangle& a, const KDTNodeTriangle& b) -> bool { return (a.tricenter[0] > b.tricenter[0]); });
		} else if (depth % 3 == 1) {
			std::sort(begin, end, [](const KDTNodeTriangle& a, const KDTNodeTriangle& b) -> bool { return (a.tricenter[1] > b.tricenter[1]); });
		} else {
			std::sort(begin, end, [](const KDTNodeTriangle& a, const KDTNodeTriangle& b) -> bool { return (a.tricenter[2] > b.tricenter[2]); });
		}

		const int midInd = (endInd - beginInd) / 2 + beginInd;
		auto mid	 = NodeList.begin();
		std::advance(mid, midInd);

		if (beginInd != midInd) {
			NodeList[midInd].rightchild		     = ConstructKDTree(beginInd, midInd, depth + 1);
			NodeList[NodeList[midInd].rightchild].parent = midInd;

			if (depth % 3 == 0) {
				NodeList[midInd].realboundright = std::min_element(begin, mid, [](const auto& a, const auto& b) -> bool { return (a.AABBmin[0] < b.AABBmin[0]); })->AABBmin[0];
			} else if (depth % 3 == 1) {
				NodeList[midInd].realboundright = std::min_element(begin, mid, [](const auto& a, const auto& b) -> bool { return (a.AABBmin[1] < b.AABBmin[1]); })->AABBmin[1];
			} else {
				NodeList[midInd].realboundright = std::min_element(begin, mid, [](const auto& a, const auto& b) -> bool { return (a.AABBmin[2] < b.AABBmin[2]); })->AABBmin[2];
			}
		}

		if ((midInd + 1) != endInd) {
			NodeList[midInd].leftchild		    = ConstructKDTree(midInd + 1, endInd, depth + 1);
			NodeList[NodeList[midInd].leftchild].parent = midInd;

			auto midp = mid + 1;
			if (depth % 3 == 0) {
				NodeList[midInd].realboundleft = std::max_element(midp, end, [](const auto& a, const auto& b) -> bool { return (a.AABBmax[0] < b.AABBmax[0]); })->AABBmax[0];
			} else if (depth % 3 == 1) {
				NodeList[midInd].realboundleft = std::max_element(midp, end, [](const auto& a, const auto& b) -> bool { return (a.AABBmax[1] < b.AABBmax[1]); })->AABBmax[1];
			} else {
				NodeList[midInd].realboundleft = std::max_element(midp, end, [](const auto& a, const auto& b) -> bool { return (a.AABBmax[2] < b.AABBmax[2]); })->AABBmax[2];
			}
		}

		return midInd;
	}

	std::vector<unsigned int> SearchNearest(const VEC Query)
	{
		double tempmax = std::numeric_limits<SCALAR>::max();
		double tempmin = std::numeric_limits<SCALAR>::max();

		std::vector<iss> ClosestCandidates;

		SearchNearestRecursive(root, Query, ClosestCandidates, tempmax, tempmin, 0);

		std::vector<unsigned int> ret;

		for (unsigned int i = 0; i < ClosestCandidates.size(); i++)
			ret.push_back(NodeList[ClosestCandidates[i].index].TriInd);

		return ret;
	}
};

}

#endif
