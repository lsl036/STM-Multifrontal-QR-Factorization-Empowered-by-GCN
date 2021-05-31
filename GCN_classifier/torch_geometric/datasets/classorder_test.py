import json
import os.path as osp

import torch
import numpy as np
import scipy.sparse as sp
from google_drive_downloader import GoogleDriveDownloader as gdd
from torch_geometric.data import Dataset,InMemoryDataset, Data


class ClassorderTest(InMemoryDataset):
    r"""The Classorder dataset from the `"GraphSAINT: Graph Sampling Based
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
        super(ClassorderTest, self).__init__(root, transform, pre_transform)
        path = self.processed_paths[0]
        if train_type==0:
            path = self.processed_paths[0]  
        
        self.data, self.slices = torch.load(path)
    @property
    def raw_file_names(self):
        return ['classorder.cites', 'classorder.content']

    @property
    def processed_file_names(self):
        return ['{}.pt'.format(s) for s in ['test']]

    def download(self):
        pass 

    def process(self):
        # 从文件中读取数据
        fcites = np.loadtxt(osp.join(self.raw_dir, 'classorder.cites'))
        fcontent = np.loadtxt(osp.join(self.raw_dir, 'classorder.content'))
        fylabel = np.loadtxt(osp.join(self.raw_dir, 'graph_y.txt'))
        
        # 人为划分划分数据集
        index_min = np.min(fcontent[:,0])
        data_length = int(np.max(fcontent[:,0])- index_min + 1)
        print(data_length)
        data_no = np.arange(0, data_length, 1)
        
        test_dt = data_no[0:data_length]

        test_data_list = []

        for i in range(data_length):
            # 每次选择一张图进行处理
            # edge_index_temp = fcites[:,1:2]
            edge_index_temp = fcites[np.where(fcites[:,0]-index_min == i)]
            row_tmp = edge_index_temp[:,1]
            col_tmp = edge_index_temp[:,2]
            row = torch.from_numpy(row_tmp).to(torch.long)
            col = torch.from_numpy(col_tmp).to(torch.long)
            edge_index = torch.stack([row, col], dim=0)
            
            # 每次选择一张图进行处理
            content_ori_temp = fcontent[np.where(fcontent[:,0]-index_min == i)]
            content_ori_temp = content_ori_temp[content_ori_temp[:,1].argsort()]
            # 所有数据归一化，除以最大的值
            # content_temp = content_ori_temp / content_ori_temp.max(axis=0)
            content_temp = content_ori_temp
            # print 取所有信息 入度数，出度数。
            x_temp = content_temp[:,2:5]
            # print 取所有信息 图序号，行号。
            no_temp = content_temp[:,0:2]
            # x数据归一化，除以最大的值
            x_temp = x_temp / (x_temp.max(axis=0)+1)
            # x_temp = np.expand_dims(x_temp, axis=1)
            # embeding需要long
            x_np = np.concatenate((no_temp,x_temp),axis=1)
            x = torch.from_numpy(x_np).to(torch.float)
            # x = torch.from_numpy(x_temp).to(torch.long)

            # x数据归一化，除以最大的值
            # y_temp = (y_temp+1) / (y_temp.max(axis=0)+2)
            # y_temp = 1/(y_temp+1)
            # y_temp = np.expand_dims(y_temp, axis=1)

            y_temp = fylabel[np.where(fylabel[:,0]-index_min == i)]
            y = y_temp[:,1]
            if 0 == i:
                print(y)
            y = torch.from_numpy(y).to(torch.long)
            data = Data(x=x, edge_index=edge_index, y=y)
            data = data if self.pre_transform is None else self.pre_transform(data)

            if i in test_dt:
                test_data_list.append(data)

        data = data if self.pre_transform is None else self.pre_transform(data)
        print(test_data_list)
        print(self.processed_paths[0])

        torch.save((self.collate(test_data_list)), self.processed_paths[0])
        # torch.save(self.collate([data]), self.processed_paths[0])

    def __repr__(self):
        return '{}()'.format(self.__class__.__name__)


