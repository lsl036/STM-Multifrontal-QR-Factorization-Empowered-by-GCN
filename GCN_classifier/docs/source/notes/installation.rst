Installation
============

We have outsourced a lot of functionality of PyTorch Geometric to other packages, which needs to be additionally installed.
These packages come with their own CPU and GPU kernel implementations based on `C++/CUDA extensions <https://github.com/pytorch/extension-cpp/>`_.

.. note::
    We do not recommend installation as root user on your system python.
    Please setup an `Anaconda/Miniconda <https://conda.io/docs/user-guide/install/index.html/>`_ environment or create a `Docker image <https://www.docker.com/>`_.

Please follow the steps below for a successful installation.

Installation via Binaries
-------------------------

We provide pip wheels for all major OS/PyTorch/CUDA combinations:

#. Ensure that at least PyTorch 1.4.0 is installed:

    .. code-block:: none

        $ python -c "import torch; print(torch.__version__)"
        >>> 1.7.0

#. Find the CUDA version PyTorch was installed with:

    .. code-block:: none

        $ python -c "import torch; print(torch.version.cuda)"
        >>> 10.2

#. Install the relevant packages:

    .. code-block:: none

         pip install torch-scatter -f https://pytorch-geometric.com/whl/torch-${TORCH}+${CUDA}.html
         pip install torch-sparse -f https://pytorch-geometric.com/whl/torch-${TORCH}+${CUDA}.html
         pip install torch-cluster -f https://pytorch-geometric.com/whl/torch-${TORCH}+${CUDA}.html
         pip install torch-spline-conv -f https://pytorch-geometric.com/whl/torch-${TORCH}+${CUDA}.html
         pip install torch-geometric

    where :obj:`${CUDA}` and :obj:`${TORCH}` should be replaced by your specific CUDA version (:obj:`cpu`, :obj:`cu92`, :obj:`cu101`, :obj:`cu102`, :obj:`cu110`) and PyTorch version (:obj:`1.4.0`, :obj:`1.5.0`, :obj:`1.6.0`, :obj:`1.7.0`), respectively.
    For example, for PyTorch 1.7.0/1.7.1 and CUDA 11.0, type:

    .. code-block:: none

         pip install torch-scatter -f https://pytorch-geometric.com/whl/torch-1.7.0+cu110.html
         pip install torch-sparse -f https://pytorch-geometric.com/whl/torch-1.7.0+cu110.html
         pip install torch-cluster -f https://pytorch-geometric.com/whl/torch-1.7.0+cu110.html
         pip install torch-spline-conv -f https://pytorch-geometric.com/whl/torch-1.7.0+cu110.html
         pip install torch-geometric

    For PyTorch 1.6.0 and CUDA 10.2, type:

    .. code-block:: none

         pip install torch-scatter -f https://pytorch-geometric.com/whl/torch-1.6.0+cu102.html
         pip install torch-sparse -f https://pytorch-geometric.com/whl/torch-1.6.0+cu102.html
         pip install torch-cluster -f https://pytorch-geometric.com/whl/torch-1.6.0+cu102.html
         pip install torch-spline-conv -f https://pytorch-geometric.com/whl/torch-1.6.0+cu102.html
         pip install torch-geometric

Installation from Source
------------------------

In case a specific version is not supported by `our wheels <https://pytorch-geometric.com/whl/>`_, you can alternatively install PyTorch Geometric from source:

#. Ensure that your CUDA is setup correctly (optional):

    #. Check if PyTorch is installed with CUDA support:

        .. code-block:: none

            $ python -c "import torch; print(torch.cuda.is_available())"
            >>> True

    #. Add CUDA to :obj:`$PATH` and :obj:`$CPATH` (note that your actual CUDA path may vary from :obj:`/usr/local/cuda`):

        .. code-block:: none

            $ export PATH=/usr/local/cuda/bin:$PATH
            $ echo $PATH
            >>> /usr/local/cuda/bin:...

            $ export CPATH=/usr/local/cuda/include:$CPATH
            $ echo $CPATH
            >>> /usr/local/cuda/include:...

    #. Add CUDA to :obj:`$LD_LIBRARY_PATH` on Linux and to :obj:`$DYLD_LIBRARY_PATH` on macOS (note that your actual CUDA path may vary from :obj:`/usr/local/cuda`):

        .. code-block:: none

            $ export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH
            $ echo $LD_LIBRARY_PATH
            >>> /usr/local/cuda/lib64:...

            $ export DYLD_LIBRARY_PATH=/usr/local/cuda/lib:$DYLD_LIBRARY_PATH
            $ echo $DYLD_LIBRARY_PATH
            >>> /usr/local/cuda/lib:...

    #. Verify that :obj:`nvcc` is accessible from terminal:

        .. code-block:: none

            $ nvcc --version
            >>> 10.2

    #. Ensure that PyTorch and system CUDA versions match:

        .. code-block:: none

            $ python -c "import torch; print(torch.version.cuda)"
            >>> 10.2

            $ nvcc --version
            >>> 10.2

#. Install all needed packages:

    .. code-block:: none

      pip install torch-scatter
      pip install torch-sparse
      pip install torch-cluster
      pip install torch-spline-conv
      pip install torch-geometric


In rare cases, CUDA or Python path problems can prevent a successful installation.
:obj:`pip` may even signal a successful installation, but runtime errors complain about missing modules, *.e.g.*, :obj:`No module named 'torch_*.*_cuda'`, or execution simply crashes with :obj:`Segmentation fault (core dumped)`.
We collected a lot of common installation errors in the `Frequently Asked Questions <https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html#frequently-asked-questions>`_ subsection.
In case the FAQ does not help you in solving your problem, please create an `issue <https://github.com/rusty1s/pytorch_geometric/issues>`_.
You should additionally verify that your CUDA is set up correctly by following the official `installation guide <https://docs.nvidia.com/cuda/index.html>`_, and that the `official extension example <https://github.com/pytorch/extension-cpp>`_ runs on your machine.

Frequently Asked Questions
--------------------------

#. ``ImportError: ***: cannot open shared object file: No such file or directory``: Add CUDA to your ``$LD_LIBRARY_PATH`` (see `Issue#43 <https://github.com/rusty1s/pytorch_geometric/issues/43>`_).

#. ``undefined symbol:``, *e.g.* ``_ZN2at6detail20DynamicCUDAInterface10set_deviceE``: Clear the pip cache and reinstall the respective package (see `Issue#7 <https://github.com/rusty1s/pytorch_scatter/issues/7>`_). On macOS, it may help to install clang compilers via conda (see `Issue#18 <https://github.com/rusty1s/pytorch_geometric/issues/18>`_):

   .. code-block:: none

      $ conda install -y clang_osx-64 clangxx_osx-64 gfortran_osx-64

#. Unable to import ``*_cuda``: You need to ``import torch`` first before importing any of the extension packages (see `Issue#6 <https://github.com/rusty1s/pytorch_scatter/issues/6>`_).

#. ``error: command '/usr/bin/nvcc' failed with exit status 2``: Ensure that at least CUDA >= 8 is installed (see `Issue#25a <https://github.com/rusty1s/pytorch_geometric/issues/25>`_ and `Issue#106 <https://github.com/rusty1s/pytorch_geometric/issues/106>`_).

#. ``return __and_<is_constructible<_Elements, _UElements&&>...>::value``: Ensure that your ``gcc`` version is at least 4.9 (and below 6) (see `Issue#25b <https://github.com/rusty1s/pytorch_scatter/issues/25>`_).
   You will also need to reinstall PyTorch because ``gcc`` versions must be consistent across all PyTorch packages.

#. ``file not recognized: file format not recognized``: Clean the repository and temporarily rename Anaconda's ``ld`` linker (see `Issue#16683 <https://github.com/pytorch/pytorch/issues/16683>`_).

#. ``undefined symbol: __cudaPopCallConfiguration``: Ensure that your PyTorch CUDA version and system CUDA version match (see `Issue#19 <https://github.com/rusty1s/pytorch_scatter/issues/19>`_):

   .. code-block:: none

      $ python -c "import torch; print(torch.version.cuda)"
      $ nvcc --version

#. ``undefined symbol: _ZN3c105ErrorC1ENS_14SourceLocationERKSs``: The ``std::string`` abi does not match between building PyTorch and its extensions.
   This is fixable by building extensions with ``-D_GLIBCXX_USE_CXX11_ABI=1`` or building PyTorch from source (see `this PyTorch thread <https://discuss.pytorch.org/t/undefined-symbol-when-import-lltm-cpp-extension/32627>`_).

#. On macOS: ``'gcc' failed with exit status 1``: Install the respective packages by using the following environment variables (see `Issue#21 <https://github.com/rusty1s/pytorch_scatter/issues/21>`_):

   .. code-block:: none

       $ MACOSX_DEPLOYMENT_TARGET=10.9 CC=clang CXX=clang++ python setup.py install

#. On macOS: ``ld: warning: directory not found for option '-L/usr/local/cuda/lib64'`` and ``ld: library not found for -lcudart``: Symlink ``cuda/lib`` to ``cuda/lib64`` (see `Issue#116 <https://github.com/rusty1s/pytorch_geometric/issues/116>`_):

   .. code-block:: none

       $ sudo ln -s /usr/local/cuda/lib /usr/local/cuda/lib64

#. On macOS: ``The version of the host compiler ('Apple clang') is not supported``: Downgrade your command line tools (see `this StackOverflow thread <https://stackoverflow.com/questions/36250949/revert-apple-clang-version-for-nvcc/46574116>`_) with the respective version annotated in the `CUDA Installation Guide for Mac <https://developer.download.nvidia.com/compute/cuda/10.1/Prod/docs/sidebar/CUDA_Installation_Guide_Mac.pdf>`_ (Section 1.1) for your specific CUDA version.
   You can download previous command line tool versions `here <https://idmsa.apple.com/IDMSWebAuth/signin?appIdKey=891bd3417a7776362562d2197f89480a8547b108fd934911bcbea0110d07f757&path=%2Fdownload%2Fmore%2F&rv=1>`_.
