import os.path as osp

import torch
import torch.nn.functional as F
# from torch.nn import Linear
from torch_geometric.data import DataLoader
from torch_geometric.datasets import Reorder 
# from torch_geometric.data import Reorder 
from torch_geometric import transforms as T
import numpy as np
import torch
from torch.nn import Sequential as Seq, Linear as Lin, ReLU
from torch_geometric.nn import GCNConv, TopKPooling
from torch_geometric.nn import global_add_pool
from torch_geometric.nn import TopKPooling
from torch_geometric.nn import global_mean_pool as gap, global_max_pool as gmp
from torch_scatter import scatter_mean

import time

test_dt = [1,2,3,4,7,8,11,12,13,14,15,18,20]
names_list=['LFAT5','LF10','Trefethen_20b','ex5','mesh1e1','Trefethen_300','bcsstk20','494_bus','nos7','1138_bus','bcsstk21','s3rmt3m3','msc01440','Muu','fv1','bcsstk25','Dubcova1','olafu','minsurfo','finan512','apache1']
embed_dim = 10
batch_size = 1
path = osp.join(osp.dirname(osp.realpath(__file__)), '.', 'data', 'Reorder')
dataset = Reorder(path,train_type=2)
test_loader = DataLoader(dataset=dataset, batch_size=batch_size)



class Net(torch.nn.Module):
    def __init__(self):
        super(Net, self).__init__()

        self.conv1 = GCNConv(embed_dim, 128)
        # self.pool1 = TopKPooling(128, ratio=0.8)
        self.conv2 = GCNConv(128, 128)
        # self.pool2 = TopKPooling(128, ratio=0.8)
        # self.conv3 = SAGEConv(128, 128)
        # self.pool3 = TopKPooling(128, ratio=0.8)
        self.item_embedding = torch.nn.Embedding(num_embeddings=100, embedding_dim=embed_dim)
        # self.lin1 = torch.nn.Linear(128, 128)
        # self.lin2 = torch.nn.Linear(128, 64)
        self.lin3 = torch.nn.Linear(128, 1)
        # self.bn1 = torch.nn.BatchNorm1d(128)
        # self.bn2 = torch.nn.BatchNorm1d(64)
        # self.act1 = torch.nn.ReLU()
        # self.act2 = torch.nn.ReLU()
  
    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        # print(x.size(),edge_index.size())
        x = self.item_embedding(x)
        x = x.squeeze()
        x = F.relu(self.conv1(x, edge_index))
        # print("x:",x.size())

        # x, edge_index, _, batch, _, _ = self.pool1(x, edge_index, None, batch)
        # print("px:",x.size())
        # x1 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
        # print("cx:",x.size())
        x = F.dropout(x, p=0.5, training=self.training)

        x = F.relu(self.conv2(x, edge_index))
     
        # x, edge_index, _, batch, _, _ = self.pool2(x, edge_index, None, batch)
        # x2 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        # x = F.relu(self.conv3(x, edge_index))

        # x, edge_index, _, batch, _, _ = self.pool3(x, edge_index, None, batch)
        # x3 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
        # print("cx:",x3.size())

        # x = x1 + x2 + x3
        # print("x:",x.size())

        # x = self.lin1(x)
        # x = self.act1(x)
        # x = self.lin2(x)
        # x = self.act2(x)
        x = self.lin3(x).view(-1)
        x = torch.sigmoid(x)
        return x


device = torch.device('cpu')
model = torch.load('./models/model_best_sign_ep450_0.0652.torch', map_location=lambda storage, loc: storage)
# model = Net().to(device)
print(model)
optimizer = torch.optim.Adam(model.parameters(), lr=0.05)

@torch.no_grad()
def test(loader):
    model.eval()
    # total_correct = total_examples = 0
    for data in loader:
        data = data.to(device)
        out = model(data)
    # return total_correct
       
@torch.no_grad()
def saveall(loader):
    sava_list = []
    model.eval()
    total_correct = total_examples = 0
    for i,data in enumerate(loader):
        # print(data)
        data = data.to(device)
        out = model(data)
        # print(out)
        _ , out_indices = torch.sort(out,dim=0)
        _ , y_indices = torch.sort(data.y,dim=0)
        total_correct += int((out_indices==y_indices).sum())
        order_list = list(out_indices.cpu().numpy().flatten())
        sava_list = [0] * len(order_list)
        for item,val in enumerate(order_list):
            sava_list[val] = item
        # print(i)
        # print(test_dt[i])
        # print(names_list[test_dt[i]])
        # sava_list = [i ]
        np.savetxt("./resultdata/{}.mtx".format(names_list[test_dt[i]]), sava_list, fmt='%d', delimiter='\t', newline='\n')
        # sava_list.extend(list(out_indices.cpu().numpy().flatten()))
    # print(sava_list)
    # np.savetxt("s3rmt3m3.txt", sava_list, fmt='%d', delimiter='\t', newline='\n')
    return total_correct

time_start = time.time()
test(test_loader)
time_end=time.time()
# print('totally cost',time_end-time_start)
test_acc = saveall(test_loader)
print('totally cost',time_end-time_start)




















