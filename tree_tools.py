import numpy as np
import pandas as pd
import copy


class Node:
    def __init__(self, name):
        self.name = name
        self.children = []
        self.children_names = []
        self.parent = None
        self.visited = False
        self.exp = None

    def update_children(self, child):
        """
        Add a node to children list
        :param child: a Node object
        :return: None
        """
        if child.name not in self.children_names:
            self.children.append(child)
            self.children_names.append(child.name)

    def update_parent(self, parent):
        """
        Add a parent node
        :param parent: a Node object
        :return: None
        """
        self.parent = parent

    def get_subtree_names(self):
        """
        Go over all children and create a list of node names
        :return: a list of names
        """
        names = [self.name]
        for child in self.children:
            names.extend(child.get_subtree_names())
        return names

    def get_ancestors(self):
        """
        Go over the ancestors of the current node and a create a list of node names
        :return: a list of names (including current node)
        """
        names = [self.name]
        if self.parent is not None:
            names.extend(self.parent.get_ancestors())
        return names

    def get_child(self, name):
        """
        This function receives a node name and returns the matching node in the subtree
        :param name: str (node's name)
        :return: Node
        """
        if self.name == name:
            return self
        for child in self.children:
            child_node = child.get_child(name)
            if child_node:
                return child_node
        return None

    def delete_node(self):
        """
        "Deletes" a node by updating it's children's parent to be it's parent,
        and it's parent's children's to include it's children.
        :return: None
        """
        for child in self.children:
            child.update_parent(self.parent)
            if self.parent:
                self.parent.update_children(child)
        if self.parent:
            self.parent.remove_child(self.name)

    def remove_child(self, child_name):
        """
        Removes a child from children lists
        :param child_name: str
        :return: None
        """
        child_node = self.get_child(child_name)
        if child_node in self.children:
            self.children.remove(child_node)
            self.children_names.remove(child_name)

    def is_split(self):
        """
        Checks if a node is a split point, meaning it has more than one child.
        :return: bool
        """
        if len(self.children) > 1:
            return True
        return False

    def get_non_splits_in_subtree(self):
        """
        A recursive function that returns all the non-splits (leafs and inner nodes with only one child)
        :return: a list of Nodes
        """
        if not self.is_split():
            non_splits = [self]
        else:
            non_splits = []
        for child in self.children:
            non_splits += child.get_non_splits_in_subtree()
        return non_splits

    def check_cycle(self):
        """
        A helper function used for debugging
        :return: bool
        """
        if self.visited:
            return True
        self.visited = True
        for child in self.children:
            if child.check_cycle():
                return True
        return False

    def update_name(self, new_name):
        """
        A helper function that updates the nodes' name (also in it's parent list)
        :param new_name: str
        :return: None
        """
        old_name = self.name
        self.name = new_name
        if self.parent:
            self.parent.remove_child(old_name)
            self.parent.update_children(self)

    def permute_tree(self, permute_dict):
        """
        Permute tree based on premute_dict (replacing the old name with the matching one in the dict).
        Note that we only change the node's name, but children_names list is unchanged (because the permutations are
        performed after the tree is built, when children_names list is no longer used.
        :param permute_dict: dict, keys are old names and values are new names
        :return: None
        """
        self.name = permute_dict[self.name]
        for child in self.children:
            child.permute_tree(permute_dict)

    def remove_suffix(self):
        """
        A helper function that removes the suffix (of length 1) from all the nodes in the subtree
        :return: None
        """
        old_name = self.name
        self.update_name(old_name[:-1])
        for child in self.children:
            child.remove_suffix()

    def calc_dist_to_child(self, child_name, cur_dist=0):
        """
        A recursive function that calculates the distance between the node and one of it's children
        :param child_name: str
        :param cur_dist: current dist from the real parent we're interested in.
        :return: int, distance from current node to child node + cur_dist
        """
        if self.name == child_name:
            return cur_dist
        best_dist = np.inf
        for child in self.children:
            best_dist = min(best_dist, child.calc_dist_to_child(child_name, cur_dist+1))
        return best_dist

    def find_lca(self, node1_name, node2_name):
        """
        A recursive function that finds the lowest common ancestor of 2 nodes.
        :param node1_name: first node's name
        :param node2_name: second node's name
        :return: Node object
        """
        if (not self.get_child(node1_name)) or (not self.get_child(node2_name)):
            return None
        lca = self
        # Check if one of the node's children is also a parent of both nodes, and therefore it should be returned.
        for child in self.children:
            child_lca = child.find_lca(node1_name, node2_name)
            if child_lca:
                lca = child_lca
        return lca

    def add_expression_profile(self, cell_profile_df):
        """
        A recursive function that adds expression profile to all the nodes
        :param cell_profile_df: a data_handeling frame with columns for each cell type (node) and row for each gene
        :return: None
        """
        self.exp = cell_profile_df[self.name]
        for child in self.children:
            child.add_expression_profile(cell_profile_df)

    def get_tree_df(self, df=None):
        """
        A helper function that creates a dataframe of the tree's structure
        :param df: a dataframe object, all the nodes in the subtree under the current node will be added to df
        :return: df
        """
        # If this is the root, create df
        if df is None:
            df = pd.DataFrame(data=None, columns=['name', 'parent', 'children'])

        # Created node dict
        cur_dict = {'name': self.name}
        if self.parent:
            cur_dict['parent'] = self.parent.name
        else:
            cur_dict['parent'] = None
        children = [child.name for child in self.children]
        cur_dict['children'] = ";".join(children)

        # Add current node's data_handeling to the dataframe
        cur_df = pd.DataFrame(data=cur_dict, index=[0])
        df = df.append(cur_df, ignore_index=True)

        # Add the children to the dataframe by recursively calling the function
        for child in self.children:
            df = child.get_tree_df(df)
        return df


def generate_tree(tree_csv_path, root):
    """
    A helper function to generate a Node object for the root (from a list describing the tree's structure)
    :param tree_csv_path: a csv file with tree's structure
    :param root: str, a node name to return
    :return: Node object (matching the root), dict (with all the nodes)
    """
    tree_df = pd.read_csv(tree_csv_path, index_col=0)
    nodes_dict = {}

    # Go over the tree structure file and generate the tree accordingly
    for row in tree_df.iterrows():
        row = row[1]
        node = nodes_dict.get(row['name'].upper(), Node(row['name'].upper()))
        # Update the parent to include the node in it's children, and the node's parent to be the parent node
        if not pd.isnull(row['parent']):
            parent_name = row['parent'].upper()
            parent = nodes_dict.setdefault(parent_name, Node(parent_name))
            node.update_parent(parent)
            parent.update_children(node)

        # Update all the children
        if not pd.isnull(row['children']):
            children_name = row['children'].upper()
            for child_name in children_name.split(";"):
                child = nodes_dict.setdefault(child_name, Node(child_name))
                child.update_parent(node)
                node.update_children(child)
        nodes_dict[row['name'].upper()] = node

    # Check that root matches a node, and if it has a parent - delete it
    root_node = nodes_dict.get(root.upper(), None)
    if root_node:
        root_node.update_parent(None)
    else:
        raise Exception("Root doesn't exist!")

    return root_node, nodes_dict


def generate_trimmed_tree(tree_csv_path, root, trim_factor):
    """
    A helper function that generates a full tree and then trims it according to trim factor
    :param tree_csv_path: a csv file with tree's structure
    :param root: str, a node name to return
    :param trim_factor: float, the fraction of nodes to trim (out of non splits)
    :return: Node object (matching the root)
    """
    main_tree, nodes_dict = generate_tree(tree_csv_path, root)

    if main_tree.check_cycle():
        raise Exception("Tree has cycles")
    non_splits = [node for node in nodes_dict.values() if not node.is_split()]
    trim_num = round(trim_factor * len(non_splits))
    nodes_to_remove = np.random.choice(non_splits, trim_num, replace=False)
    for node in nodes_to_remove:
        node_name = node.name
        node.delete_node()
        if node_name not in nodes_dict:
            print("a")
            continue
        nodes_dict.pop(node_name)

    return main_tree, nodes_dict


def permute_tree(root):
    """
    A helper function that permute a tree
    :param root: Node object
    :return: Node object
    """
    # Create a dict of old to new names
    root_copy = copy.deepcopy(root)

    nodes = root.get_subtree_names()
    permute_dict = {i: j for i, j in zip(nodes, np.random.permutation(nodes))}
    root_copy.permute_tree(permute_dict)
    return root
