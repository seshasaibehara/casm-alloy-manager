import setuptools

setuptools.setup(
    name="casm-alloy-manager",
    version="0.0.1",
    packages=[
        "casmam",
        "casmam.xtal",
        "casmam.scripts",
        "casmam.mapping",
        "casmam.xtallib",
        "casmam.xtallib.top5",
    ],
    package_data={"": ["*.vasp"]},
    entry_points={
        "console_scripts": ["casm-alloy-manager=casmam.scripts.casm_alloy_manager:main"]
    },
    # install_requires=["numpy", "scipy"],
    python_requires=">=3.7",
)
