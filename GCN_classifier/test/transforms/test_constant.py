import torch
from torch_geometric.transforms import Constant
from torch_geometric.data import Data


def test_constant():
    assert Constant().__repr__() == 'Constant(value=1)'

    x = torch.tensor([[-1, 0], [0, 0], [2, 0]], dtype=torch.float)
    edge_index = torch.tensor([[0, 1], [1, 2]])

    data = Data(edge_index=edge_index, num_nodes=3)
    data = Constant()(data)
    assert len(data) == 2
    assert data.edge_index.tolist() == edge_index.tolist()
    assert data.x.tolist() == [[1], [1], [1]]

    data = Data(edge_index=edge_index, x=x)
    data = Constant()(data)
    assert len(data) == 2
    assert data.edge_index.tolist() == edge_index.tolist()
    assert data.x.tolist() == [[-1, 0, 1], [0, 0, 1], [2, 0, 1]]
