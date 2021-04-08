import sys
import random
import os.path as osp
import shutil

import torch
from torchvision.datasets.mnist import MNIST, read_image_file, read_label_file
import torchvision.transforms as T
from torch_geometric.data import download_url, extract_gz, DataLoader
from torch_geometric.data.makedirs import makedirs
from torch_geometric.transforms import ToSLIC

resources = [
    'http://yann.lecun.com/exdb/mnist/t10k-images-idx3-ubyte.gz',
    'http://yann.lecun.com/exdb/mnist/t10k-labels-idx1-ubyte.gz',
]


def test_to_superpixels():
    root = osp.join('/', 'tmp', str(random.randrange(sys.maxsize)))

    raw_folder = osp.join(root, 'MNIST', 'raw')
    processed_folder = osp.join(root, 'MNIST', 'processed')

    makedirs(raw_folder)
    makedirs(processed_folder)
    for resource in resources:
        path = download_url(resource, raw_folder)
        extract_gz(path, osp.join(root, raw_folder))

    test_set = (
        read_image_file(osp.join(raw_folder, 't10k-images-idx3-ubyte')),
        read_label_file(osp.join(raw_folder, 't10k-labels-idx1-ubyte')),
    )

    torch.save(test_set, osp.join(processed_folder, 'training.pt'))
    torch.save(test_set, osp.join(processed_folder, 'test.pt'))

    dataset = MNIST(root, download=False)

    dataset.transform = T.Compose([T.ToTensor(), ToSLIC()])

    data, y = dataset[0]
    assert len(data) == 2
    assert data.pos.dim() == 2 and data.pos.size(1) == 2
    assert data.x.dim() == 2 and data.x.size(1) == 1
    assert data.pos.size(0) == data.x.size(0)
    assert y == 7

    loader = DataLoader(dataset, batch_size=2, shuffle=False)
    for data, y in loader:
        assert len(data) == 4
        assert data.pos.dim() == 2 and data.pos.size(1) == 2
        assert data.x.dim() == 2 and data.x.size(1) == 1
        assert data.batch.dim() == 1
        assert data.ptr.dim() == 1
        assert data.pos.size(0) == data.x.size(0) == data.batch.size(0)
        assert y.tolist() == [7, 2]
        break

    dataset.transform = T.Compose(
        [T.ToTensor(), ToSLIC(add_seg=True, add_img=True)])

    data, y = dataset[0]
    assert len(data) == 4
    assert data.pos.dim() == 2 and data.pos.size(1) == 2
    assert data.x.dim() == 2 and data.x.size(1) == 1
    assert data.pos.size(0) == data.x.size(0)
    assert data.seg.size() == (1, 28, 28)
    assert data.img.size() == (1, 1, 28, 28)
    assert data.seg.max().item() + 1 == data.x.size(0)
    assert y == 7

    loader = DataLoader(dataset, batch_size=2, shuffle=False)
    for data, y in loader:
        assert len(data) == 6
        assert data.pos.dim() == 2 and data.pos.size(1) == 2
        assert data.x.dim() == 2 and data.x.size(1) == 1
        assert data.batch.dim() == 1
        assert data.ptr.dim() == 1
        assert data.pos.size(0) == data.x.size(0) == data.batch.size(0)
        assert data.seg.size() == (2, 28, 28)
        assert data.img.size() == (2, 1, 28, 28)
        assert y.tolist() == [7, 2]
        break

    shutil.rmtree(root)
