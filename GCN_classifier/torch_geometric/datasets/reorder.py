import json
import os.path as osp

import torch
import numpy as np
import scipy.sparse as sp
from google_drive_downloader import GoogleDriveDownloader as gdd
from torch_geometric.data import Dataset,InMemoryDataset, Data


class Reorder(InMemoryDataset):
    r"""The Reorder dataset from the `"GraphSAINT: Graph Sampling Based
    Inductive Learning Method" <https://arxiv.org/abs/1907.04931>`_ paper,
    containing descriptions and common properties of images.

    Args:
        root (string): Root directory where the dataset should be saved.
        transform (callable, optional): A function/transform that takes in an
            :obj:`torch_geometric.data.Data` object and returns a transformed
            version. The data object will be transformed before every access.
            (default: :obj:`None`)
        pre_transform (callable, optional): A function/transform that takes in
            an :obj:`torch_geometric.data.Data` object and returns a
            transformed version. The data object will be transformed before
            being saved to disk. (default: :obj:`None`)
    """

    # adj_full_id = '1crmsTbd1-2sEXsGwa2IKnIB7Zd3TmUsy'
    # feats_id = '1join-XdvX3anJU_MLVtick7MgeAQiWIZ'
    # class_map_id = '1uxIkbtg5drHTsKt-PAsZZ4_yJmgFmle9'
    # role_id = '1htXCtuktuCW8TR8KiKfrFDAxUgekQoV7'

    def __init__(self, root,
                 train_type=0, transform=None, pre_transform=None):
        super(Reorder, self).__init__(root, transform, pre_transform)
        path = self.processed_paths[0]
        if train_type==0:
            path = self.processed_paths[0]  
        elif train_type==1:
            path= self.processed_paths[1]
        elif train_type==2:
            path= self.processed_paths[2]
        self.data, self.slices = torch.load(path)
    @property
    def raw_file_names(self):
        return ['reorder.cites', 'reorder.content']

    @property
    def processed_file_names(self):
        return ['{}.pt'.format(s) for s in ['train', 'val', 'test']]

    def download(self):
        pass 

    def process(self):
        # 从文件中读取数据
        fcites = np.loadtxt(osp.join(self.raw_dir, 'reorder.cites'))
        fcontent = np.loadtxt(osp.join(self.raw_dir, 'reorder.content'))
        # train_dt = [5,6,10]
        # val_dt = [8,9]
        # test_dt = [0,1,2,3,4,7,11,12,13,14,15]
        train_dt = [5,6,10,16,17,19]
        val_dt = [0,9]
        test_dt = [1,2,3,4,7,8,11,12,13,14,15,18,20]
        train_data_list, val_data_list, test_data_list = [], [], []
        for i in range(20):
            # 每次选择一张图进行处理
            # edge_index_temp = fcites[:,1:2]
            edge_index_temp = fcites[np.where(fcites[:,0] == i)]
            row_tmp = edge_index_temp[:,1]
            col_tmp = edge_index_temp[:,2]
            row = np.append(row_tmp, col_tmp)
            col = np.append(col_tmp, row_tmp)
            row = torch.from_numpy(row).to(torch.long)
            col = torch.from_numpy(col).to(torch.long)
            edge_index = torch.stack([row, col], dim=0)
            
            # 每次选择一张图进行处理
            content_ori_temp = fcontent[np.where(fcontent[:,0] == i)]
            content_ori_temp = content_ori_temp[content_ori_temp[:,1].argsort()]
            # 所有数据归一化，除以最大的值
            # content_temp = content_ori_temp / content_ori_temp.max(axis=0)
            content_temp = content_ori_temp
            # print 只取度的信息
            x_temp = content_temp[:,2:3]
            # x数据归一化，除以最大的值
            # x_temp = x_temp / x_temp.max(axis=0)
            # x_temp = np.expand_dims(x_temp, axis=1)
            # embeding需要long
            # x = torch.from_numpy(x_temp).to(torch.float)
            x = torch.from_numpy(x_temp).to(torch.long)

            y_temp = content_temp[:,-1]
            # x数据归一化，除以最大的值
            y_temp = (y_temp+1) / (y_temp.max(axis=0)+2)
            # if 0 == i:
            #     print(y_temp)
            # y_temp = 1/(y_temp+1)
            # y_temp = np.expand_dims(y_temp, axis=1)
            y = torch.from_numpy(y_temp).to(torch.float)
            data = Data(x=x, edge_index=edge_index, y=y)
            data = data if self.pre_transform is None else self.pre_transform(data)
            if i in train_dt:
                train_data_list.append(data)
            if i in val_dt:
                val_data_list.append(data)
            if i in test_dt:
                test_data_list.append(data)
            # role_list = np.arange(0, x_temp.shape[0])
            # np.random.shuffle(role_list)

            # role_tr = role_list[0 : int(x_temp.shape[0]*0.2)]
            # role_val = role_list[int(x_temp.shape[0]*0.2) : int(x_temp.shape[0]*0.5)]
            # role_test = role_list[0 : int(x_temp.shape[0])]

            # train_mask = torch.zeros(x.size(0), dtype=torch.bool)
            # train_mask[torch.tensor(role_tr)] = True

            # val_mask = torch.zeros(x.size(0), dtype=torch.bool)
            # val_mask[torch.tensor(role_val)] = True

            # test_mask = torch.zeros(x.size(0), dtype=torch.bool)
            # test_mask[torch.tensor(role_test)] = True

        # data = Data(x=x, edge_index=edge_index, y=y, train_mask=train_mask,
        #             val_mask=val_mask, test_mask=test_mask)

        data = data if self.pre_transform is None else self.pre_transform(data)
        torch.save((self.collate(train_data_list)), self.processed_paths[0])
        torch.save((self.collate(val_data_list)), self.processed_paths[1])
        torch.save((self.collate(test_data_list)), self.processed_paths[2])
        # torch.save(self.collate([data]), self.processed_paths[0])

    def __repr__(self):
        return '{}()'.format(self.__class__.__name__)


    # def get(self, index):
    #     # print(self.slices['edge_index'])
    #     x0, x1 = self.slices['x'][index], self.slices['x'][index+1]
    #     y0, y1 = self.slices['y'][index], self.slices['y'][index+1]
    #     edge_index0, edge_index1 = self.slices['edge_index'][index], self.slices['edge_index'][index+1]
    #     x = self.data['x'][x0: x1]
    #     y = self.data['y'][y0: y1]
    #     edge_index = self.data['edge_index'][:,edge_index0: edge_index1]
    #     # data = torch.load(osp.join(self.processed_dir, 'data_{}.pt'.format(index)))
    #     data = Data(x=x, edge_index=edge_index, y=y)
    #     return data

    # def len(self):
    #     # print("len:",self.slices['x'].size(0))
    #     return self.slices['x'].size(0)-1

    # def process(self):
    #     # 从文件中读取数据
    #     # f = open(osp.join(self.raw_dir, 'reorder.cites'), 'r', encoding="utf-8")
    #     # for lines in f.readlines():
    #     #     line_data = lines.split()
    #     # f.close()
    #     fcites = np.loadtxt(osp.join(self.raw_dir, 'reorder.cites'))
    #     fcontent = np.loadtxt(osp.join(self.raw_dir, 'reorder.content'))

    #     # f = np.load(osp.join(self.raw_dir, 'adj_full.npz'))
    #     # adj = sp.csr_matrix((f['data'], f['indices'], f['indptr']), f['shape'])
    #     # adj = adj.tocoo()
    #     # row = torch.from_numpy(adj.row).to(torch.long)
    #     # col = torch.from_numpy(adj.col).to(torch.long)
    #     # 暂时只选择一张图
    #     # edge_index_temp = fcites[:,1:2]
    #     edge_index_temp = fcites[np.where(fcites[:,0] == 11)]
    #     # edge_index_temp = edge_index_temp[:,1:3]
    #     row_tmp = edge_index_temp[:,1]
    #     col_tmp = edge_index_temp[:,2]
    #     row = np.append(row_tmp, col_tmp)
    #     col = np.append(col_tmp, row_tmp)
    #     row = torch.from_numpy(row).to(torch.long)
    #     col = torch.from_numpy(col).to(torch.long)
    #     # print(row, col)
    #     edge_index = torch.stack([row, col], dim=0)
    #     # edge_index = edge_index_1
    #     # print(edge_index_1)
    #     # edge_index_2 = torch.stack([col, row], dim=0)
    #     # edge_index = np.concatenate(edge_index_1, edge_index_2)
    #     # # edge_index_temp[:, [0,1]] = edge_index_temp[:, [1,0]]
    #     # print(edge_index[0],edge_index[1])
    #     # edge_index_temp = fcites[np.where(fcites[:,0] == 0)]
    #     # edge_index = torch.from_numpy(edge_index_temp).to(torch.long)
        
    #     # 暂时只选择一张图
    #     content_ori_temp = fcontent[np.where(fcontent[:,0] == 11)]
    #     content_temp = content_ori_temp / content_ori_temp.max(axis=0)
    #     # print
    #     x_temp = content_temp[:,0:-1]
    #     # x_temp = x_temp[np.where(x_temp[:,0] == 0)]
    #     # print(x.shape())
    #     x = torch.from_numpy(x_temp).to(torch.float)

    #     # content_temp = fcontent[np.where(fcontent[:,0] == 0)]
    #     y_temp = content_temp[:,-1]
    #     # y_max = max(y_temp)
    #     y_temp = np.expand_dims(y_temp, axis=1)
    #     # ys = np.zeros((int(y_max+1),int(y_max+1)))
    #     # for i,v in enumerate(y_temp):
    #     #     ys[i,int(v)] = 1
    #     # y = torch.from_numpy(ys).to(torch.long)
    #     y = torch.from_numpy(y_temp).to(torch.float)

    #     # ys = [-1] * x.size(0)
    #     # with open(osp.join(self.raw_dir, 'class_map.json')) as f:
    #     #     class_map = json.load(f)
    #     #     for key, item in class_map.items():
    #     #         ys[int(key)] = item
    #     # y = torch.tensor(ys)

    #     # with open(osp.join(self.raw_dir, 'role.json')) as f:
    #     #     role = json.load(f)
    #     role_list = np.arange(0, x_temp.shape[0])
    #     np.random.shuffle(role_list)

    #     role_tr = role_list[0 : int(x_temp.shape[0]*0.2)]
    #     role_val = role_list[int(x_temp.shape[0]*0.2) : int(x_temp.shape[0]*0.5)]
    #     # role_test = role_list[int(x_temp.shape[0]*0.5) : int(x_temp.shape[0])]
    #     role_test = role_list[0 : int(x_temp.shape[0])]
    #     print(len(role_test))

    #     # role_tr = np.random.randint(0,x.size(0)-1,(int(x.size(0)*0.1)))
    #     # role_val = np.random.randint(0,x.size(0)-1,(int(x.size(0)*0.3)))
    #     # role_test = np.random.randint(0,x.size(0)-1,(int(x.size(0)*0.6)))

    #     train_mask = torch.zeros(x.size(0), dtype=torch.bool)
    #     train_mask[torch.tensor(role_tr)] = True

    #     val_mask = torch.zeros(x.size(0), dtype=torch.bool)
    #     val_mask[torch.tensor(role_val)] = True

    #     test_mask = torch.zeros(x.size(0), dtype=torch.bool)
    #     test_mask[torch.tensor(role_test)] = True

    #     data = Data(x=x, edge_index=edge_index, y=y, train_mask=train_mask,
    #                 val_mask=val_mask, test_mask=test_mask)

    #     data = data if self.pre_transform is None else self.pre_transform(data)

    #     torch.save(self.collate([data]), self.processed_paths[0])
