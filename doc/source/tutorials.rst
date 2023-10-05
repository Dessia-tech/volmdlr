Tutorials
=========

In this example, a regular decision tree with next_unique_node is implemented
to find the bast solution for a Knapsack Problem. A selection of items among a
collection of 10 items, having a mass and a value, have to be made to yield
the highest value while fitting in a knapsack of maximum 30kg. To solve this
problem the regular decision tree will be composed of 10 levels (a selection of
maximum 10 items) and each node will have 10 possibilities, one for each item.
However an item can't be picked twice so this example uses next_unique_node. ::

  from dectree import RegularDecisionTree


  # Creating a class for the items
  class Item:
      def __init__(self, mass, value, name):
          self.mass, self.value, self.name = mass, value, name

      def __repr__(self):
          return f'{self.name} : value = {self.value} ; mass = {self.mass}'


  # Knapsack maximum weight
  max_mass = 30

  # List of all the items
  items = [Item(4, 10, 'item1'),
           Item(5, 11, 'item2'),
           Item(6, 14, 'item3'),
           Item(7, 18, 'item4'),
           Item(8, 20, 'item5'),
           Item(9, 24, 'item6'),
           Item(10, 27, 'item7'),
           Item(11, 28, 'item8'),
           Item(12, 30, 'item9'),
           Item(13, 33, 'item10')]

  # Creating a regular decision tree
  tree = RegularDecisionTree(np=[len(items)] * len(items))

  result = None
  max_value = 0
  while not tree.finished:
      # Stocking information about current node
      current_mass = sum(items[i].mass for i in tree.current_node)
      current_value = sum(items[i].value for i in tree.current_node)
      current_items = [items[i] for i in tree.current_node]

      # Checking viability of current node
      valid = False if current_mass > max_mass else True

      if valid and current_value > max_value:
          result = current_items
          max_value = current_value

      # Going to next node
      tree.next_unique_node(current_node_viability=valid)

However, knowing that order of the items inside a solution doesn't matter,
the algorithm can be improved using next_sorted_unique_node instead of
next_unique_node. It is then equivalent to a regular decision tree with
10 levels, on for each item, and only 2 possibilities for each node : take
the corresponding item or not. ::

  # Creating a regular decision tree
  tree = RegularDecisionTree(np=[2] * len(items))

  result = None
  max_value = 0
  while not tree.finished:
      # Stocking information about current node
      current_mass = sum(items[item_i].mass * i for item_i, i in enumerate(tree.current_node))
      current_value = sum(items[item_i].value * i for item_i, i in enumerate(tree.current_node))
      current_items = [items[item_i] for item_i, i in enumerate(tree.current_node) if i]

      # Checking viability of current node
      valid = False if current_mass > max_mass else True

      if valid and current_value > max_value:
          result = current_items
          max_value = current_value

      # Going to next node
      tree.next_node(current_node_viability=valid)
