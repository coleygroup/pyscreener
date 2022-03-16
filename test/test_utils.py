import numpy as np
import pytest

from pyscreener.utils import Reduction, reduce_scores


@pytest.fixture(params=[(10,), (50,), (10, 1), (100, 5), (1000, 10)])
def size(request):
    return request.param


@pytest.fixture
def s(size):
    return np.random.normal(size=size)


@pytest.fixture(params=[1, 5, 10])
def k(request):
    return request.param


@pytest.fixture(params=list(Reduction))
def reduction(request):
    return request.param


@pytest.fixture(params=[1, "foo", "best", 3.14])
def invalid_reduction(request):
    return request.param


@pytest.mark.parametrize("invalid_reduction,", [(1,), ("foo",), ("best",), (3.14,)])
def test_reduce_invalid_reduction(invalid_reduction):
    s = np.empty(10)
    with pytest.raises(ValueError):
        reduce_scores(s, invalid_reduction)


def test_reduce_scores_best(s):
    np.testing.assert_array_equal(reduce_scores(s, Reduction.BEST), s.min(-1))


def test_reduce_scores_avg(s):
    np.testing.assert_array_equal(reduce_scores(s, Reduction.AVG), s.mean(-1))


def test_reduce_scores_boltz(s):
    s_actual = reduce_scores(s, Reduction.BOLTZMANN)
    s_e = np.exp(-s)
    Z = np.sum(s_e, -1, keepdims=True)
    s_desired = np.sum(s * s_e / Z, -1)

    np.testing.assert_allclose(s_actual, s_desired, rtol=1e-9, atol=1e-6)


def test_reduce_scores_topk(s, k):
    s_actual = reduce_scores(s, Reduction.TOP_K, k=k)
    s_desired = np.sort(s, -1)[..., :k].mean(-1)

    np.testing.assert_array_equal(s_actual, s_desired)


def test_reduce_all_nan(size, reduction):
    s = np.empty(size)
    s[:] = np.nan

    np.testing.assert_array_equal(reduce_scores(s, reduction), s.sum(-1))
