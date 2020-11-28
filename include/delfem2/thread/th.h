/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <functional>
#include <mutex>
#include <queue>
#include <thread>

#ifndef DFM2_TH_H
#define DFM2_TH_H

namespace delfem2 {
namespace thread {

unsigned int mymin(unsigned int x, unsigned int y) {
  return (x <= y) ? x : y;
}

unsigned int mymax(unsigned int x, unsigned int y) {
  return (x >= y) ? x : y;
}

template<typename FUNCTION>
void parallel_for(
    unsigned int ntask,
    FUNCTION function,
    unsigned int target_concurrency = 0) {
  if( ntask == 0 ){ return; }
  unsigned int nthread = (target_concurrency == 0) ? std::thread::hardware_concurrency()
                                                               : target_concurrency;
  nthread = mymin(ntask, (nthread == 0) ? 4 : nthread);
  const unsigned int ntasks_per_thread = (ntask / nthread) + (ntask % nthread == 0 ? 0 : 1);
  auto tasks_for_each_thread = [&](const unsigned int ithread) {
    const unsigned int itask_start = ithread * ntasks_per_thread;
    const unsigned int itask_end = mymin(ntask, itask_start + ntasks_per_thread);
    for (unsigned int itask = itask_start; itask < itask_end; ++itask) {
      function(itask);
    }
  };
  std::vector<std::thread> aThread;
  for (unsigned int jthread = 0; jthread < nthread; ++jthread) {
    aThread.push_back(std::thread(tasks_for_each_thread, jthread));
  }
  for (auto &t : aThread) { t.join(); }
}

}
}

#endif /* TH_H */