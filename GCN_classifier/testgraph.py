import os.path as osp

import torch
import torch.nn.functional as F
# from torch.nn import Linear
from torch_geometric.data import DataLoader
from torch_geometric.datasets import ClassorderTest
from torch_geometric.datasets import Classorder 
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


batch_size = 1


path = osp.join(osp.dirname(osp.realpath(__file__)), '.', 'data', 'Classtest')
# path = osp.join(osp.dirname(osp.realpath(__file__)), '.', 'data', 'Classorder')
dataset = ClassorderTest(path,train_type=0)
# dataset = Classorder(path,train_type=2)
test_loader = DataLoader(dataset=dataset, batch_size=batch_size)

# 额外信息文件读取并保存 xelist=[]
xelist = np.loadtxt('./data/Classtest/raw/QR_extinfo.txt')
# xelist = np.loadtxt('./data/Classorder/raw/QR_extinfo.txt')
def getexinfo(graph_no_list):
    result = None
    for graph_no in graph_no_list:
        result_temp = xelist[np.where(xelist[:,0] == graph_no)]
        if result is None:
            result = result_temp[:,1:]
        else:
            # result = result + result_temp
            result = np.concatenate((result,result_temp[:,1:]),axis=0)
        # result.append(result_temp[:,1:].astype(np.float32)) 

    
    # print(result)
    result = torch.from_numpy(np.array(result).astype(np.float32))
    return result.to(device)

class Net(torch.nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        # 卷积神经网络
        self.conv1 = GraphConv(train_dataset.num_features-2, 128)
        self.pool1 = TopKPooling(128, ratio=0.8)
        self.conv2 = GraphConv(128, 128)
        self.pool2 = TopKPooling(128, ratio=0.8)

        # self.conv3 = GraphConv(128, 128)
        # self.pool3 = TopKPooling(128, ratio=0.8)

        # 加上2遍conv再加上额外信息
        self.lin1 = torch.nn.Linear(128*2+10, 64)
        self.lin2 = torch.nn.Linear(64, 32)
        self.lin3 = torch.nn.Linear(32, train_dataset.num_classes)

    def forward(self, data, exinfo):
        x, edge_index, batch = data.x[:,2:5], data.edge_index, data.batch
        x1 = F.relu(self.conv1(x, edge_index))
        x1 = F.dropout(x1, p=0.2, training=self.training)
        x2 = F.relu(self.conv2(x1, edge_index))
        x2 = F.dropout(x2, p=0.2, training=self.training)

        x = torch.cat([x1, x2], dim=-1)
        x = gmp(x, batch)
        xe = exinfo
        x = torch.cat([x, xe], dim=-1)

        x = F.relu(self.lin1(x))
        x = F.dropout(x, p=0.2, training=self.training)
        x = F.relu(self.lin2(x))
        x = F.log_softmax(self.lin3(x), dim=-1)

        return x

# device = torch.device('cpu')
# model = torch.load('./models/model_nodecluster_ep270_0.0812.torch', map_location=lambda storage, loc: storage)

device = torch.device('cuda:0')
model = torch.load('./models/graphmodel_ep5000.torch')

# model = Net().to(device)
print(model)
# optimizer = torch.optim.Adam(model.parameters(), lr=0.05)

@torch.no_grad()
def test(loader):
    model.eval()
    total_correct = total_examples = 0
    
    for data in loader:
        graph_id = data.x[:,0].cpu().numpy()
        exinfo = getexinfo(list(set(data.x[:,0].numpy())))
        data = data.to(device)
        
        for val_t in data.x[:,0].cpu().numpy():
            id_print = []
            id_print.append(val_t)
        time_start = time.time()
        out = model(data, exinfo)
        time_end = time.time()

        pred = out.max(1)[1]
        print("id=", id_print, "out=", pred, "y=", data.y, ", time=", time_end - time_start)
        total_correct += int((pred==data.y).sum())
        total_examples += int((pred==pred).sum())
    

    return total_correct, total_examples


total_correct, total_examples = test(test_loader)
print('正确数%d,总数%d: '%(total_correct,total_examples))
print('accuracy: ',total_correct/total_examples)



















