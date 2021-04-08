import os.path as osp

import torch
import torch.nn.functional as F
# from torch.nn import Linear
from torch_geometric.data import DataLoader
from torch_geometric.datasets import Classorder 
# from torch_geometric.data import Classorder 
from torch_geometric import transforms as T
import numpy as np
import torch
from torch.nn import Sequential as Seq, Linear as Lin, ReLU
from torch_geometric.nn import GraphConv, TopKPooling
from torch_geometric.nn import global_add_pool
from torch_geometric.nn import TopKPooling
from torch_geometric.nn import global_mean_pool as gap, global_max_pool as gmp
from torch_scatter import scatter_mean

import time

embed_dim = 10
batch_size = 1
num_epochs = 3000

path = osp.join(osp.dirname(osp.realpath(__file__)), '.', 'data', 'Classorder')
dataset = Classorder(path,train_type=2)
test_loader = DataLoader(dataset=dataset, batch_size=batch_size)


class Net(torch.nn.Module):
    def __init__(self):
        super(Net, self).__init__()

        self.conv1 = GraphConv(train_dataset.num_features, 128)
        self.pool1 = TopKPooling(128, ratio=0.8)
        self.conv2 = GraphConv(128, 128)
        self.pool2 = TopKPooling(128, ratio=0.8)
        # self.conv3 = GraphConv(128, 128)
        # self.pool3 = TopKPooling(128, ratio=0.8)

        self.lin1 = torch.nn.Linear(256, 128)
        self.lin2 = torch.nn.Linear(128, 64)
        self.lin3 = torch.nn.Linear(64, train_dataset.num_classes)

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch

        x = F.relu(self.conv1(x, edge_index))
        x, edge_index, _, batch, _, _ = self.pool1(x, edge_index, None, batch)
        x1 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        x = F.relu(self.conv2(x, edge_index))
        x, edge_index, _, batch, _, _ = self.pool2(x, edge_index, None, batch)
        x2 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        # x = F.relu(self.conv3(x, edge_index))
        # x, edge_index, _, batch, _, _ = self.pool3(x, edge_index, None, batch)
        # x3 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        # x = x1 + x2 + x3
        x = x1 + x2

        x = F.relu(self.lin1(x))
        x = F.dropout(x, p=0.5, training=self.training)
        x = F.relu(self.lin2(x))
        x = F.log_softmax(self.lin3(x), dim=-1)

        return x

# device = torch.device('cpu')
# model = torch.load('./models/model_nodecluster_ep270_0.0812.torch', map_location=lambda storage, loc: storage)

device = torch.device('cuda:0')
model = torch.load('./models/graphmodel_ep50.torch')

# model = Net().to(device)
print(model)
# optimizer = torch.optim.Adam(model.parameters(), lr=0.05)

@torch.no_grad()
def test(loader):
    model.eval()
    total_correct = total_examples = 0
    for data in loader:
        data = data.to(device)
        out = model(data)
        pred = out.max(1)[1]
        # print("out",pred,"y",data.y)
        total_correct += int((pred==data.y).sum())
        total_examples += int((pred==pred).sum())
        
    return total_correct, total_examples
        # total_examples += idx.numel()

# @torch.no_grad()
# def saveall(loader):
#     sava_list = []
#     model.eval()
#     total_correct = total_examples = 0
#     # print(loader[1])
#     for i,data in enumerate(loader):
#         # print(data.x.numpy())
#         # print(data)
#         data = data.to(device)
#         out = model(data)
#         # print(out)
#         _ , out_indices = torch.sort(out,dim=0)
#         # _ , y_indices = torch.sort(data.y[:,0],dim=0)
#         # total_correct += int((out_indices==y_indices).sum())
#         order_list = list(out_indices.cpu().numpy().flatten())
#         # sava_list = [0] * len(order_list)
#         # for item,val in enumerate(order_list):
#         #     sava_list[val] = item
#         # print(i)
#         # print(test_dt[i])
#         # print(names_list[test_dt[i]])
#         # sava_list = [i]
#         filename = int(data.x.cpu().numpy()[0,0])
#         np.savetxt("./resultdata/{}.mtx".format(filename), order_list, fmt='%d', delimiter='\t', newline='\n')
#         # sava_list.extend(list(out_indices.cpu().numpy().flatten()))
#     # print(sava_list)
#     # np.savetxt("s3rmt3m3.txt", sava_list, fmt='%d', delimiter='\t', newline='\n')
#     return total_correct

time_start = time.time()
total_correct, total_examples = test(test_loader)
# test_acc = saveall(test_loader)
time_end=time.time()
# print('totally cost',time_end-time_start)
print('totally cost: ',time_end-time_start)
print('正确数%d,总数%d: '%(total_correct,total_examples))
print('accuracy: ',total_correct/total_examples)



















