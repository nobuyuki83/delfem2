/*
 * Copyright (c) 2020 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef DFM2_THREAD_H
#define DFM2_THREAD_H

#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <thread>

#include "delfem2/dfm2_inline.h"

namespace delfem2 {

template<typename T, typename Func>
inline void parallel_for(
    T num,
    Func &&func,
    [[maybe_unused]] unsigned int target_concurrency = 0) {
  auto futures = std::vector<std::future<void>>{};

  const unsigned int nthreads = (target_concurrency == 0) ?
                                std::thread::hardware_concurrency() : target_concurrency;
  std::atomic<T> next_idx(0);
  std::atomic<bool> has_error(false);
  for (auto thread_id = 0; thread_id < (int) nthreads; thread_id++) {
    futures.emplace_back(
        std::async(std::launch::async, [&func, &next_idx, &has_error, num]() {
          try {
            while (true) {
              auto idx = next_idx.fetch_add(1);
              if (idx >= num) break;
              if (has_error) break;
              func(idx);
            }
          } catch (...) {
            has_error = true;
            throw;
          }
        }));
  }
  for (auto &f: futures) f.get();
}

template<typename T, typename Func>
inline void parallel_for(
    T num1,
    T num2,
    Func &&func,
    unsigned int target_concurrency = 0) {
  auto futures = std::vector<std::future<void>>{};
  const unsigned int nthreads = (target_concurrency == 0) ?
                                std::thread::hardware_concurrency() : target_concurrency;
  std::atomic<T> next_idx(0);
  std::atomic<bool> has_error(false);
  for (auto thread_id = 0; thread_id < (int) nthreads; thread_id++) {
    futures.emplace_back(std::async(
        std::launch::async, [&func, &next_idx, &has_error, num1, num2]() {
          try {
            while (true) {
              auto j = next_idx.fetch_add(1);
              if (j >= num2) break;
              if (has_error) break;
              for (auto i = (T) 0; i < num1; i++) func(i, j);
            }
          } catch (...) {
            has_error = true;
            throw;
          }
        }));
  }
  for (auto &f: futures) f.get();
}

}

#endif /* DFM2_THREAD_H */