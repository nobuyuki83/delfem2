####################################################################
# Copyright (c) 2019 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################


import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--is_gl", action="store", default="true", help="my option: type1 or type2"
    )

'''
@pytest.fixture
def cmdopt(request):
    return request.config.getoption("--is_gl")
'''