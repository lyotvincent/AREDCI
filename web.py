from flask import  Flask, render_template, request
import sys, os, argparse

ROOT_PATH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT_PATH, 'src'))
from src.main import run_dci


app=Flask(__name__)
app.config['DEBUG'] = True  # 开启 debug

@app.route('/')
def index():
    msg = "Welcome to AREDCI"
    return render_template("index.html", data=msg)

@app.route('/run_dci', methods=['POST'])
def do_dci():
    # 从请求中获取参数
    files = request.json['files']
    samples = request.json['samples']
    groups = request.json['groups']
    gap = int(request.json['gap'])
    blacklist_path = request.json['blacklist_path']
    remove_loops_in_blacklist = request.json['remove_loops_in_blacklist']
    remove_self_ligated_loops = request.json['remove_self_ligated_loops']
    fdr_threshold = float(request.json['fdr_threshold'])
    loop_len_threshold = int(request.json['loop_len_threshold'])
    intra_only = request.json['intra_only']
    chr_filter = request.json['chr_filter']
    pet_threshold = int(request.json['pet_threshold'])
    norm_method = request.json['norm_method']
    fig_demand = request.json['fig_demand']
    reproducibility_checkbox = request.json['reproducibility_checkbox']
    idr_checkbox = request.json['idr_checkbox']
    dci_method = request.json['dci_method']
    output_path = request.json['output_path']
    print("files:", files)
    print("samples:", samples)
    print("groups:", groups)
    print("gap:", gap)
    print("blacklist_path:", blacklist_path)
    print("remove_loops_in_blacklist:", remove_loops_in_blacklist)
    print("remove_self_ligated_loops:", remove_self_ligated_loops)
    print("fdr_threshold:", fdr_threshold)
    print("intra_only:", intra_only)
    print("chr_filter:", type(chr_filter))
    print("chr_filter:", chr_filter)
    print("pet_threshold:", pet_threshold)
    print("norm_method:", norm_method)
    print("fig_demand:", fig_demand)
    print("reproducibility_checkbox:", reproducibility_checkbox)
    print("idr_checkbox:", idr_checkbox)
    print("dci_method:", dci_method)
    print("output_path:", output_path)

    reproducibility_result = run_dci(files=files,
            samples=samples,
            groups=groups,
            blacklist_path=blacklist_path,
            gap=gap,
            remove_loops_in_blacklist=remove_loops_in_blacklist,
            remove_self_ligated_loops=remove_self_ligated_loops,
            fdr_threshold=fdr_threshold,
            intra_only=intra_only,
            loop_len_threshold=loop_len_threshold,
            chr_filter=chr_filter,
            pet_threshold=pet_threshold,
            norm_method=norm_method,
            fig_demand=fig_demand,
            reproducibility_checkbox=reproducibility_checkbox,
            idr_checkbox=idr_checkbox,
            dci_method=dci_method,
            output_path=output_path)

    # 返回响应
    return {'success': True, 'data': 'temp none', 'aredci_result': reproducibility_result}

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Run the web server")
    parser.add_argument('-p', '--port', type=int, default=1080, help='The port to run the server on')
    args = parser.parse_args()

    app.run(port=args.port, host="127.0.0.1")

