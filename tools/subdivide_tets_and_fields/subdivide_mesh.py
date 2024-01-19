import subprocess

def run_shell_script(mesh_path, fra_path, subd_mesh_path):
    print(12)
    # Path to the shell script
    script_path = './call_subdivide_mesh.sh'

    args = [script_path, mesh_path, fra_path, subd_mesh_path]

    print(args)

    # Call the script with subprocess.Popen
    try:
        result = subprocess.Popen(args,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = result.communicate()
        if result.returncode == 0:
            print("Script executed successfully.")
            print("Output:", output.decode())
        else:
            print("Script execution failed.")
            print("Error:", error.decode())
    except Exception as e:
        print("An error occurred:", e)


run_shell_script("~/Downloads/spot_5000_sing.mesh", "~/Downloads/frames_curr.fra", "~/Downloads/outtest.mesh")

# Example usage
# call_subdivide_mesh('/path/to/mesh', '/path/to/fra', '/path/to/subd_mesh')