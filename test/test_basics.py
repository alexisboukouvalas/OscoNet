'''
Basic test installation is working
'''
import collections

import numpy as np
import pytest

from OscopeBootstrap import SyntheticDataset
from OscopeBootstrap.oscope_tf import bootstrap_hypothesis_test

DataParameters = collections.namedtuple('DataParameters', 'ngroups NG N G noiselevel')


@pytest.fixture
def dataset_options() -> DataParameters:
    return DataParameters(
        ngroups=3,
        NG=3,  # number of genes in group
        G=100,  # total number of genes
        N=5,  # number of cells
        noiselevel=1,
    )

@pytest.fixture
def dataset_generation(dataset_options):
    data, phaseG, angularSpeed = SyntheticDataset.GetSimISyntheticData(NG=dataset_options.NG,
                                                                       G=dataset_options.G,
                                                                       N=dataset_options.N,
                                                                       noiseLevel=dataset_options.noiselevel,
                                                                       ngroups=dataset_options.ngroups)
    return data, phaseG, angularSpeed


def test_generation(dataset_options, dataset_generation):
    data, phaseG, angularSpeed = dataset_generation
    assert data.shape == (dataset_options.G, dataset_options.N)
    assert phaseG.shape == (dataset_options.G, )  # phase for each gene
    assert angularSpeed.shape == (dataset_options.G, )
    assert np.all(~np.isnan(data))
