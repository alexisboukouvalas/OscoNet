from OscopeBootstrap.SyntheticDataset import GetSimISyntheticData


def test_data_generation_single_group():
    data, _, _ = GetSimISyntheticData(NG=2, G=10, N=100, noiseLevel=3, ngroups=1)
    G, N = data.shape
    assert G == 10
    assert N == 100


def test_data_generation_three_groups():
    data, _, _ = GetSimISyntheticData(NG=5, G=100, N=100, noiseLevel=3, ngroups=3)
    G, N = data.shape
    assert G == 100
    assert N == 100
