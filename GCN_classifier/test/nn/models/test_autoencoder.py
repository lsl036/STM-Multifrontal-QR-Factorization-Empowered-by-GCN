import torch

from torch import Tensor as T
from torch_geometric.nn import GAE, VGAE, ARGA, ARGVA
from torch_geometric.data import Data
from torch_geometric.utils import train_test_split_edges


def test_gae():
    model = GAE(encoder=lambda x: x)
    model.reset_parameters()

    x = torch.Tensor([[1, -1], [1, 2], [2, 1]])
    z = model.encode(x)
    assert z.tolist() == x.tolist()

    adj = model.decoder.forward_all(z)
    assert adj.tolist() == torch.sigmoid(
        torch.Tensor([[+2, -1, +1], [-1, +5, +4], [+1, +4, +5]])).tolist()

    edge_index = torch.tensor([[0, 1], [1, 2]])
    value = model.decode(z, edge_index)
    assert value.tolist() == torch.sigmoid(torch.Tensor([-1, 4])).tolist()

    edge_index = torch.tensor([[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                               [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]])
    data = Data(edge_index=edge_index)
    data.num_nodes = edge_index.max().item() + 1
    data = train_test_split_edges(data, val_ratio=0.2, test_ratio=0.3)

    z = torch.randn(11, 16)
    loss = model.recon_loss(z, data.train_pos_edge_index)
    assert loss.item() > 0

    auc, ap = model.test(z, data.val_pos_edge_index, data.val_neg_edge_index)
    assert auc >= 0 and auc <= 1 and ap >= 0 and ap <= 1


def test_vgae():
    model = VGAE(encoder=lambda x: (x, x))

    x = torch.Tensor([[1, -1], [1, 2], [2, 1]])
    model.encode(x)
    assert model.kl_loss().item() > 0

    model.eval()
    model.encode(x)


def test_arga():
    model = ARGA(encoder=lambda x: x, discriminator=lambda x: T([0.5]))
    model.reset_parameters()

    x = torch.Tensor([[1, -1], [1, 2], [2, 1]])
    z = model.encode(x)

    assert model.reg_loss(z).item() > 0
    assert model.discriminator_loss(z).item() > 0


def test_argva():
    model = ARGVA(encoder=lambda x: (x, x), discriminator=lambda x: T([0.5]))

    x = torch.Tensor([[1, -1], [1, 2], [2, 1]])
    model.encode(x)
    model.reparametrize(model.__mu__, model.__logstd__)
    assert model.kl_loss().item() > 0


def test_init():
    encoder = torch.nn.Linear(16, 32)
    decoder = torch.nn.Linear(32, 16)
    discriminator = torch.nn.Linear(32, 1)

    GAE(encoder, decoder)
    VGAE(encoder, decoder)
    ARGA(encoder, discriminator, decoder)
    ARGVA(encoder, discriminator, decoder)
