import json

with open("input_dir/standardized/242698.json", "r") as json_file:
    sample_dict = json.load(json_file)

ranges = sample_dict["ranges"]
query = ""

for pos in sample_dict.keys():
    if pos == "ranges":
        continue
    query += f"{pos}{sample_dict[pos]['alt'].replace('gap', 'del')} "

print()
print(ranges)
print()
print(query)
