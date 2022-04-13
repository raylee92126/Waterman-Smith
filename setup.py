import setuptools

setuptools.setup(
        name="ws-alignment",
        version="0.0.1",
        author="Raymond Lee",
        author_email="raymond92126@gmail.com",
        description="Run Waterman-Smith Alignment Algorithm",
        long_description="This is an assignment for S&DS3532, where we coded the Waterman-Smith Algorithm with the added feature of supporting affine gaps.",
        long_description_content_type="text/markdown",
        url="https://github.com/raymond92126/Waterman-Smith",
        license="MIT",
        packages=["ws-alignment"],
        install_requires=["numpy", "pandas", "argparse"]
    )
