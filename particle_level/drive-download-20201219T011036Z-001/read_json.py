import json

with open('nbkit_example_power.json') as json_data:
    d = json.load(json_data)


print(d.keys())
print(d['edges'])
