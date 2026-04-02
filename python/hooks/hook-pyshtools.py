from PyInstaller.utils.hooks import collect_data_files, collect_submodules
import os

datas = collect_data_files('pyshtools')
hiddenimports = collect_submodules('pyshtools')

# Brute force the entire doc directory
pyshtools_path = os.path.dirname(__import__('pyshtools').__file__)
doc_path = os.path.join(pyshtools_path, 'doc')
if os.path.exists(doc_path):
    datas.append((doc_path, 'pyshtools/doc'))