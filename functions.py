import subprocess
import re
import os


def run(cmd):
    completed = subprocess.run(["powershell", "-Command", cmd], capture_output=True)
    return completed


def bb_pynomad(var):
    try:
        print(var.__dir__())
        x = [var.get_coord(i) for i in range(var.size())]  # convert the Nomad input to list

        # modify the launching script
        with open("Launch.py", "w") as f:
            f.write("from fem_model import *\n")
            f.write("from post_process import *\n")
            f.write("model(" + str(x) + ")\n")  # compute the model
            f.write("post_process(" + str(x) + ")\n")  # extract the data
            f.close()

        # launch the script to the powershell
        cmd_output = run("abaqus cae nogui=Launch.py")
        result = re.search("stdout=b'(.*)', stderr", str(cmd_output))
        print(cmd_output)
        # delete all the computation files
        for file in os.listdir(os.getcwd()):
            if file.startswith("waterbomb") or file.startswith("abaqus"):
                os.remove(os.path.join(file))

        # read the FEM result located in the Report.txt file
        energy = report(result)

        # Print the result in the python console
        print("result:" + str(-energy))

        # Send the result to Nomad
        var.setBBO(str(energy).encode("UTF-8"))
        return 1
    except Exception as e:
        print(f"An error occurred: {e}")


def report(result):
    with open("Report.txt", "r") as file:
        last_line = file.readlines()[-1]
        file.close()
        last_line = last_line[last_line.rfind("\t") + len("\t"):]
    if result.group(1) == "" and last_line != '-inf\n':
        value = -float(last_line)
    else:
        value = float('inf')
    return value
