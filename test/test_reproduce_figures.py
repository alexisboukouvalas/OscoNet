import pytest
import os
from casestudy.Whitfield.reproduce_pseudotime_figures import reproduce_figures


def get_load_dir(request, load_directory):
    file_path = os.path.join(request.fspath.dirname, load_directory)
    print(file_path)
    return file_path


@pytest.mark.parametrize("configuration", [1, 2])
def test_reproduce_figures_smoketest(request, configuration: int):
    """
    Smoketest to check we can produce both sets of figures for the different clusters
    """
    dir_path = get_load_dir(request, '../casestudy/Whitfield/Results')
    reproduce_figures(configuration, path_to_results=dir_path)
