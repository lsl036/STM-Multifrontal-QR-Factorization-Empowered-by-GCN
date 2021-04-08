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
from torch_geometric.nn import GraphConv, TopKPooling, SAGEConv
from torch_geometric.nn import global_add_pool
from torch_geometric.nn import TopKPooling
from torch_geometric.nn import global_mean_pool as gap, global_max_pool as gmp
from torch_scatter import scatter_mean

import time

embed_dim = 10
batch_size = 4
num_epochs = 10000
lr=0.0001

path = osp.join(osp.dirname(osp.realpath(__file__)), '.', 'data', 'Classorder')
train_dataset = Classorder(path,train_type=0)
val_dataset = Classorder(path,train_type=1)

# train_dataset = reorder_dataset.data
# print(train_dataset.data,train_dataset.slices)
train_loader = DataLoader(dataset=train_dataset, batch_size=batch_size)
val_loader = DataLoader(dataset=val_dataset, batch_size=batch_size)
# dataset_val = Classorder(path,train_type=1)


# val_loader = DataLoader(dataset_val[0], batch_size=32 * 128)

# data = dataset[0]
# print(data.x)
# print(dataset.num_classes)
# 额外信息文件读取并保存 xelist=[]
xelist = np.loadtxt('./data/Classorder/raw/QR_extinfo.txt')
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
        # 只使用大图的属性做MLP-NN
        # self.lin1 = torch.nn.Linear(10, 64)
        self.lin2 = torch.nn.Linear(64, 32)
        self.lin3 = torch.nn.Linear(32, train_dataset.num_classes)
# 添加额外信息8个
    def forward(self, data, exinfo):
        # 两层卷积层
        x, edge_index, batch = data.x[:,2:5], data.edge_index, data.batch
        x1 = F.relu(self.conv1(x, edge_index))
        x1 = F.dropout(x1, p=0.2, training=self.training)
        x2 = F.relu(self.conv2(x1, edge_index))
        x2 = F.dropout(x2, p=0.2, training=self.training)

        x = torch.cat([x1, x2], dim=-1)
        x = gmp(x, batch)
        xe = exinfo
        x = torch.cat([x, xe], dim=-1)

        # 直接MLP
        # x = exinfo
        x = F.relu(self.lin1(x))
        x = F.dropout(x, p=0.2, training=self.training)
        x = F.relu(self.lin2(x))
        x = F.log_softmax(self.lin3(x), dim=-1)

        return x


def graphtrain():
    model.train()

    loss_all = 0
    for data in train_loader:
        # 额外信息
        # print(list(set(data.x[:,0].numpy())))
        exinfo = getexinfo(list(set(data.x[:,0].numpy())))
        data = data.to(device)
        optimizer.zero_grad()
        output = model(data, exinfo)
        label = data.y.to(device)
        # if loss_all==0:
        #     print("out",output,"y",data.y)
        loss = F.nll_loss(output, data.y)
        # loss = crit(output, label)
        loss.backward()
        loss_all += data.num_graphs * loss.item()
        optimizer.step()
        
    return loss_all / len(train_dataset)
@torch.no_grad()
def graphtest(loader):
    model.eval()
    total_correct = total_examples = 0
    total_Traincorrect = total_Trainexamples = 0
    # time_start = time.time()
    for data in loader:
        
        # 额外信息
        exinfo = getexinfo(list(set(data.x[:,0].numpy())))
        data = data.to(device)
        out = model(data, exinfo)
        pred = out.max(1)[1]
        
        # 输出预测的信息
        # graph_id = data.x[:,0].cpu().numpy()
        # for i_t in set(data.x[:,0].cpu().numpy()):
        # 调整数据集
        # val_save = 9999
        # id_print = []
        # for val_t in data.x[:,0].cpu().numpy():
        #     if val_save != val_t:
        #         val_save = val_t
        #         id_print.append( val_save)

        # print("id=",id_print,"out=",pred,"y=",data.y)
        total_correct += int((pred==data.y).sum())
        total_examples += int((pred==pred).sum())
    # time_end=time.time()
    # print('totally cost',time_end-time_start,'s')

    # 看看训练集上的准确率
    for data_train in train_loader:
        # 额外信息
        exinfo_train = getexinfo(list(set(data_train.x[:,0].numpy())))
        data_train = data_train.to(device)
        out_train = model(data_train, exinfo_train)
        pred_train = out_train.max(1)[1]
        total_Traincorrect += int((pred_train==data_train.y).sum())
        total_Trainexamples += int((pred_train==pred_train).sum())

    return total_correct, total_examples, total_Traincorrect, total_Trainexamples
    # return total_correct, total_examples


# device = torch.device('cpu')
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
model = Net().to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=lr)

# crit = torch.nn.BCELoss()
# train_loader = DataLoader(train_dataset, batch_size=batch_size)
for epoch in range(num_epochs):
    loss = graphtrain()
    if epoch % (num_epochs/10) == 0:
        torch.save(model, './models/graphmodel_ep{}.torch'.format(epoch))
#     # train_acc = graphtest(dataset_val)
#     # val_acc = graphtest(dataset_val)
#     # # time_start = time.time()
#     # # test_acc = graphtest(test_loader)
#     # # time_end=time.time()
#     # # print('totally cost',time_end-time_start)
#     # # test_acc = saveall(test_loader)
#     # # val_acc = graphtest(val_loader)
#     # # test_acc = graphtest(test_loader)
        # train_acc = graphtest(train_loader)
    if epoch % 10 == 0:
        # total_correct, total_examples = graphtest(val_loader)
        total_correct, total_examples, total_Traincorrect, total_Trainexamples = graphtest(val_loader)
        print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, '
        # print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, Train: {train_acc:.4f}, '
                f'total_correct: {total_correct:.1f},'
                f'total_examples: {total_examples:.1f},'
                f'accuracy: {total_correct/total_examples:.4f}'
                f'Train_accuracy: {total_Traincorrect/total_Trainexamples:.4f}')




