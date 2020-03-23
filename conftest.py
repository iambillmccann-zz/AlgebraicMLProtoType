import pytest
"""
When pytest discovers a conftest.py, it modifies sys.path so it can import stuff from the conftest module. 
So, pytest will be forced to append it to sys.path. The fixture is used to setup a global test environment
for all the training_data classes.
"""
@pytest.fixture(scope = 'class')
def training_data(request):

    request.cls.test_this = 'Test this!'
    yield
