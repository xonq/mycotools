def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "integration: live network test (downloads JGI/NCBI data; needs credentials)",
    )
