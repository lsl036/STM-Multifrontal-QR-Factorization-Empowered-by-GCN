# import os.path as osp

# import torch
# import torch.nn.functional as F
# # from torch.nn import Linear
# from torch_geometric.data import DataLoader
# from torch_geometric.datasets import Reorder 
# # from torch_geometric.data import Reorder 
# from torch_geometric import transforms as T
# import numpy as np
# import torch
# from torch.nn import Sequential as Seq, Linear as Lin, ReLU
# from torch_geometric.nn import GCNConv, TopKPooling
# from torch_geometric.nn import global_add_pool
# from torch_geometric.nn import TopKPooling
# from torch_geometric.nn import global_mean_pool as gap, global_max_pool as gmp
# from torch_scatter import scatter_mean

# import time

# embed_dim = 10
# batch_size = 1
# num_epochs = 500

# path = osp.join(osp.dirname(osp.realpath(__file__)), '.', 'data', 'Reorder')
# train_dataset = Reorder(path,train_type=0)
# val_dataset = Reorder(path,train_type=1)

# # train_dataset = reorder_dataset.data
# # print(train_dataset.data,train_dataset.slices)
# train_loader = DataLoader(dataset=train_dataset, batch_size=batch_size)
# val_loader = DataLoader(dataset=val_dataset, batch_size=batch_size)
# # dataset_val = Reorder(path,train_type=1)


# # val_loader = DataLoader(dataset_val[0], batch_size=32 * 128)

# # data = dataset[0]
# # print(data.x)
# # print(dataset.num_classes)


# class Net(torch.nn.Module):
#     def __init__(self):
#         super(Net, self).__init__()

#         self.conv1 = GCNConv(embed_dim, 128)
#         # self.pool1 = TopKPooling(128, ratio=0.8)
#         self.conv2 = GCNConv(128, 128)
#         # self.pool2 = TopKPooling(128, ratio=0.8)
#         # self.conv3 = SAGEConv(128, 128)
#         # self.pool3 = TopKPooling(128, ratio=0.8)
#         self.item_embedding = torch.nn.Embedding(num_embeddings=100, embedding_dim=embed_dim)
#         # self.lin1 = torch.nn.Linear(128, 128)
#         # self.lin2 = torch.nn.Linear(128, 64)
#         self.lin3 = torch.nn.Linear(128, 1)
#         # self.bn1 = torch.nn.BatchNorm1d(128)
#         # self.bn2 = torch.nn.BatchNorm1d(64)
#         # self.act1 = torch.nn.ReLU()
#         # self.act2 = torch.nn.ReLU()
  
#     def forward(self, data):
#         x, edge_index, batch = data.x, data.edge_index, data.batch
#         x = self.item_embedding(x)
#         x = x.squeeze()        
#         # print(x,edge_index)
#         x = F.relu(self.conv1(x, edge_index))
#         # print("x:",x.size())

#         # x, edge_index, _, batch, _, _ = self.pool1(x, edge_index, None, batch)
#         # print("px:",x.size())
#         # x1 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
#         # print("cx:",x.size())
#         x = F.dropout(x, p=0.5, training=self.training)

#         x = F.relu(self.conv2(x, edge_index))
     
#         # x, edge_index, _, batch, _, _ = self.pool2(x, edge_index, None, batch)
#         # x2 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

#         # x = F.relu(self.conv3(x, edge_index))

#         # x, edge_index, _, batch, _, _ = self.pool3(x, edge_index, None, batch)
#         # x3 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
#         # print("cx:",x3.size())

#         # x = x1 + x2 + x3
#         # print("x:",x.size())

#         # x = self.lin1(x)
#         # x = self.act1(x)
#         # x = self.lin2(x)
#         # x = self.act2(x)
#         x = self.lin3(x).view(-1)
#         x = torch.sigmoid(x)
#         return x

# def train():
#     model.train()

#     loss_all = 0
#     for data in train_loader:
#         data = data.to(device)
#         optimizer.zero_grad()
#         output = model(data)
#         label = data.y.to(device)
#         # if loss_all==0:
#             # print("out",output,"y",data.y)
#         loss = F.mse_loss(output, data.y)
#         # loss = crit(output, label)
#         loss.backward()
#         loss_all += data.num_graphs * loss.item()
#         optimizer.step()
        
#     return loss_all / len(train_dataset)
# @torch.no_grad()
# def test(loader):
#     model.eval()
#     total_correct = total_examples = 0
#     for data in loader:
#         data = data.to(device)
#         out = model(data)
#         _ , out_indices = torch.sort(out,dim=0,descending=False)
#         _ , y_indices = torch.sort(data.y,dim=0,descending=False)
#         if total_correct==0:
#             print("out",out,"\ny",data.y)
#             order_list = list(out_indices.cpu().numpy().flatten())
#             sava_list = [0] * len(order_list)
#             for item,val in enumerate(order_list):
#                 sava_list[val] = item
#             y_list = list(y_indices.cpu().numpy().flatten())
#             sava_y_list = [0] * len(y_list)
#             for item,val in enumerate(y_list):
#                 sava_y_list[val] = item
#             print("out",sava_list,"\ny",sava_y_list)
#         total_correct += int((out_indices==y_indices).sum())
#     return total_correct
#         # total_examples += idx.numel()


# # device = torch.device('cpu')
# device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
# model = Net().to(device)
# optimizer = torch.optim.Adam(model.parameters(), lr=0.0002)

# # crit = torch.nn.BCELoss()
# # train_loader = DataLoader(train_dataset, batch_size=batch_size)
# for epoch in range(num_epochs):
#     loss = train()
#     if epoch % (num_epochs/10) == 0:
#         torch.save(model, './models/model_best_sign_ep{}_{:.4f}.torch'.format(epoch, loss))
# #     # train_acc = test(dataset_val)
# #     # val_acc = test(dataset_val)
# #     # # time_start = time.time()
# #     # # test_acc = test(test_loader)
# #     # # time_end=time.time()
# #     # # print('totally cost',time_end-time_start)
# #     # # test_acc = saveall(test_loader)
# #     # # val_acc = test(val_loader)
# #     # # test_acc = test(test_loader)
#         # train_acc = test(train_loader)
#     if epoch % 10 == 0:
#         val_acc = test(val_loader)
#         print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, '
#         # print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, Train: {train_acc:.4f}, '
#                 f'Val: {val_acc:.4f}')













#     # def __init__(self, in_channels=dataset.num_features):
#     #     super(Net, self).__init__()

#     #     self.conv1 = GINConv(Seq(Lin(in_channels, 64), ReLU(), Lin(64, 64)))
#     #     # self.pool1 = TopKPooling(in_channels, min_score=0.05)
#     #     self.conv2 = GINConv(Seq(Lin(64, 64), ReLU(), Lin(64, 64)))

#     #     self.lin = torch.nn.Linear(64, 1)

#     # def forward(self, data):
#     #     x, edge_index = data.x, data.edge_index
#     #     # print(x.size(), edge_index.size())
#     #     out = F.relu(self.conv1(x, edge_index))
#     #     # print(out)

#     #     # out, edge_index, _, batch, perm, score = self.pool1(
#     #     #     out, edge_index, None, batch, attn=x)
#     #     ratio = out.size(0) / x.size(0)

#     #     # out = F.relu(self.conv2(out, edge_index))
#     #     # out = global_add_pool(out, batch)
#     #     out = self.lin(out).view(-1)
#     #     return out.sigmoid()

 


# # device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# # model = Net().to(device)
# # # print(model)
# # optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

# # def train(loader):
# #     model.train()

# #     total_loss = total_examples = 0
# #     total_correct = total_examples = 0
# #     for i in range(len(loader)):
        
# #         data = loader.get(i)
# #         # train_loader = DataLoader(dataset=dataset, batch_size=4, shuffle=True)
# #         # print(data)
# #         # return 
# #     # for i,data in enumerate(loader):
# #         # dataarray = data.edge_index.numpy()
# #         # print(max(dataarray.flatten()))
# #         # print(dataarray[0])
# #         data = data.to(device)
# #         optimizer.zero_grad()
# #         out = model(data)
        
# #         # _ , out_indices = torch.sort(out,dim=0)
# #         # _ , y_indices = torch.sort(data.y,dim=0)
# #         # total_correct += int((out_indices == y_indices).sum())
# #         # print(total_correct)
# #         # print("out",i,out.cpu().detach().numpy())
# #         # print("y",data.y.cpu().numpy())
# #         loss = F.mse_loss(out, data.y)
# #         loss.backward()
# #         optimizer.step()
# #         # total_loss += float(loss) * idx.numel()
# #         # total_examples += idx.numel()

# #     # return total_loss / total_examples
# #     return float(loss)


# # @torch.no_grad()
# # def test(loader):
# #     model.eval()

# #     total_correct = total_examples = 0
# #     for i in range(len(loader)):
# #         data = loader.get(i)
# #         data = data.to(device)
# #         out = model(data)
# #         _ , out_indices = torch.sort(out,dim=0)
# #         _ , y_indices = torch.sort(data.y,dim=0)
# #         total_correct += int(((out_indices-y_indices < 10) and (out_indices-y_indices > 0)).sum())
# #     return total_correct
# #         # total_examples += idx.numel()


# # for epoch in range(1, 10):
# #     loss = train(dataset)
# #     val_acc = test(dataset_val)
# #     if epoch % 5 == 1:
# #         torch.save(model, './models/model_best_sign_ep{}_{:.4f}.torch'.format(epoch, val_acc))
# #     # train_acc = test(dataset_val)
# #     # val_acc = test(dataset_val)
# #     # # time_start = time.time()
# #     # # test_acc = test(test_loader)
# #     # # time_end=time.time()
# #     # # print('totally cost',time_end-time_start)
# #     # # test_acc = saveall(test_loader)
# #     # # val_acc = test(val_loader)
# #     # # test_acc = test(test_loader)
# #         train_acc = test(dataset)
# #         val_acc = test(dataset_val)
# #         print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, Train: {train_acc:.4f}, '
# #                 f'Val: {val_acc:.4f}')


# # def train():
# #     model.train()

# #     total_loss = total_examples = 0
# #     for idx, data in enumerate(train_loader):
# #     # for idx in train_loader:
# #         xs = [data.x[idx].to(device)]
# #         xs += [data[f'x'][idx].to(device) for i in range(1, K + 1)]
# #         y = data.y[idx].to(device)

# #         optimizer.zero_grad()
# #         out = model(xs)
# #         print(out.size(), y.size())
# #         loss = F.mse_loss(out, y)
# #         loss.backward()
# #         optimizer.step()

# #         total_loss += float(loss) * idx.numel()
# #         total_examples += idx.numel()

# #     return total_loss / total_examples


# # @torch.no_grad()
# # def test(loader):
# #     model.eval()

# #     total_correct = total_examples = 0
# #     for idx in loader:
# #         xs = [data.x[idx].to(device)]
# #         xs += [data[f'x'][idx].to(device) for i in range(1, K + 1)]
# #         y = data.y[idx].to(device)

# #         out = model(xs)
# #         # print("out",out.cpu().numpy())
# #         _ , out_indices = torch.sort(out,dim=0)
# #         _ , y_indices = torch.sort(y,dim=0)
# #         # np.savetxt("s3rmt3m3.txt", out_indices.cpu().numpy(), fmt='%d', delimiter='\t', newline='\n')
# #         # print("out_indices",out_indices.cpu().numpy())
# #         # print("y_indices",y_indices.cpu().numpy())
# #         total_correct += int((out_indices == y_indices).sum())
# #         total_examples += idx.numel()

#     # return total_correct / total_examples


