import setuptools

setuptools.setup(
    name="casm-alloy-manager",
    version="0.0.1",
    packages=[
        "casmam",
        "casmam.scripts",
        "casmam.xtallib",
        "casmam.xtallib.common",
    ],
    package_data={"": ["*.vasp"]},
    entry_points={
        "console_scripts": ["casm-alloy-manager=casmam.scripts.casm_alloy_manager:main"]
    },
    install_requires=["numpy", "pandas", "libcasm-xtal", "libcasm-mapping"],
    python_requires=">=3.7",
)
