from setuptools import setup
from setuptools_rust import Binding, RustExtension

setup(
    name="rustiq",
    version="1.0",
    author="Simon Martiel and Timothée Goubault de Brugière",
    rust_extensions=[RustExtension("rustiq.rustiq", binding=Binding.PyO3, debug=False)],
    packages=["rustiq"],
    # rust extensions are not zip safe, just like C-extensions.
    zip_safe=False,
    install_requires=[
        "numpy",
    ],
)
