#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
# from .. import rdmcl
from .. import helpers
from .. import launch_worker
import os
import sqlite3
import pandas as pd
from io import StringIO
from buddysuite import buddy_resources as br
from buddysuite import AlignBuddy as Alb
from buddysuite import SeqBuddy as Sb
from math import ceil
import time
import argparse


def mock_valueerror(*args, **kwargs):
    raise ValueError(args, kwargs)


def mock_keyboardinterupt(*args, **kwargs):
    raise KeyboardInterrupt(args, kwargs)


def test_instantiate_worker():
    temp_dir = br.TempDir()
    worker = launch_worker.Worker(temp_dir.path)
    assert worker.working_dir == temp_dir.path
    assert worker.wrkdb_path == os.path.join(temp_dir.path, "work_db.sqlite")
    assert worker.hbdb_path == os.path.join(temp_dir.path, "heartbeat_db.sqlite")
    assert worker.output == os.path.join(temp_dir.path, ".worker_output")
    assert worker.masterclear_path == os.path.join(temp_dir.path, "MasterClear")
    assert worker.heartrate == 60
    assert type(worker.heartbeat) == launch_worker.rdmcl.HeartBeat
    assert worker.heartbeat.hbdb_path == worker.hbdb_path
    assert worker.heartbeat.pulse_rate == worker.heartrate
    assert worker.heartbeat.thread_type == "worker"
    assert worker.max_wait == 120
    assert worker.cpus == br.cpu_count() - 1
    assert worker.worker_file == ""
    assert worker.data_file == ""
    assert worker.start_time
    assert worker.split_time == 0
    assert worker.idle == 1
    assert worker.running == 1
    assert worker.last_heartbeat_from_master == 0
    assert worker.subjob_num == 1
    assert worker.num_subjobs == 1
    assert worker.job_id_hash is None
    worker.heartbeat.end()


def test_start_worker_no_master(hf, capsys):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    temp_dir.copy_to("%sheartbeat_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path, heartrate=1, max_wait=1)
    with pytest.raises(SystemExit):
        worker.start()
    assert worker.split_time != 0
    assert worker.last_heartbeat_from_master != 0
    worker.data_file = ".Worker_1.dat"
    out, err = capsys.readouterr()
    assert "Starting Worker_2" in out


def test_start_worker_clean_dead_master(hf, capsys, monkeypatch):
    monkeypatch.setattr(launch_worker, "random", lambda *_: 0.991)
    # No need to actually call the function, just confirm we can get there
    monkeypatch.setattr(launch_worker.Worker, "clean_dead_threads", lambda *_: print("PASSED"))
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    temp_dir.copy_to("%sheartbeat_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path, heartrate=1, max_wait=1)
    with pytest.raises(SystemExit):
        worker.start()
    out, err = capsys.readouterr()
    assert "PASSED" in out


def test_start_worker_fetch_queue(hf, capsys):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    temp_dir.copy_to("%sheartbeat_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path, heartrate=1, max_wait=1)

    work_con = sqlite3.connect(os.path.join(temp_dir.path, "work_db.sqlite"))
    work_cursor = work_con.cursor()
    work_cursor.execute("INSERT INTO "
                        "waiting (hash, master_id) "
                        "VALUES ('foo', 2)")
    work_cursor.execute("INSERT INTO "
                        "queue (hash, psi_pred_dir, master_id, align_m, align_p, trimal, gap_open, gap_extend) "
                        "VALUES ('foo', ?, 2, 'clustalo', '', 'gappyout 50 90 clean', 0, 0)",
                        (os.path.join(hf.resource_path, "psi_pred"),))

    work_con.commit()

    seqbuddy = hf.get_data("cteno_panxs")
    seqbuddy = Sb.pull_recs(seqbuddy, "Oma")  # Only 4 records, which means 6 comparisons
    seqbuddy.write(os.path.join(worker.output, "foo.seqs"))

    with pytest.raises(SystemExit):
        worker.start()

    out, err = capsys.readouterr()
    assert "Running foo" in out
    assert "Creating MSA (4 seqs)" in out
    assert "Trimal (4 seqs)" in out
    assert os.path.isfile(os.path.join(worker.output, "foo.aln"))
    assert "Updating 4 psipred dataframes" in out
    assert "Preparing all-by-all data" in out
    assert "Running all-by-all data (6 comparisons)" in out
    assert "Processing final results" in out


def test_start_worker_1seq_error(hf):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    temp_dir.copy_to("%sheartbeat_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path, heartrate=1, max_wait=1)

    work_con = sqlite3.connect(os.path.join(temp_dir.path, "work_db.sqlite"))
    work_cursor = work_con.cursor()
    work_cursor.execute("INSERT INTO "
                        "waiting (hash, master_id) "
                        "VALUES ('foo', 2)")
    work_cursor.execute("INSERT INTO "
                        "queue (hash, psi_pred_dir, master_id, align_m, align_p, trimal, gap_open, gap_extend) "
                        "VALUES ('foo', ?, 2, 'clustalo', '', 'gappyout 50 90 clean', 0, 0)",
                        (os.path.join(hf.resource_path, "psi_pred"),))
    work_con.commit()

    # Only a single sequence present
    seqbuddy = hf.get_data("cteno_panxs")
    seqbuddy = Sb.pull_recs(seqbuddy, "Oma-PanxαC")
    seqbuddy.write(os.path.join(worker.output, "foo.seqs"))

    with pytest.raises(ValueError) as err:
        worker.start()

    assert "Queued job of size 1 encountered: foo" in str(err)


def test_start_worker_deleted_check(hf, capsys, monkeypatch):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    temp_dir.copy_to("%sheartbeat_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path, heartrate=1, max_wait=1)
    monkeypatch.setattr(helpers, "dummy_func", lambda *_: os.remove(worker.worker_file))

    with pytest.raises(SystemExit):
        worker.start()

    out, err = capsys.readouterr()
    assert "Terminating Worker_2 because of deleted check file" in out, print(out)
    assert not os.path.isfile(worker.data_file)


def test_start_worker_deal_with_subjobs(hf, capsys, monkeypatch):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    temp_dir.copy_to("%sheartbeat_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path, heartrate=1, max_wait=1)

    work_con = sqlite3.connect(os.path.join(temp_dir.path, "work_db.sqlite"))
    work_cursor = work_con.cursor()
    work_cursor.execute("INSERT INTO "
                        "waiting (hash, master_id) "
                        "VALUES ('foo', 2)")
    work_cursor.execute("INSERT INTO "
                        "queue (hash, psi_pred_dir, master_id, align_m, align_p, trimal, gap_open, gap_extend) "
                        "VALUES ('foo', ?, 2, 'clustalo', '', 'gappyout 50 90 clean', 0, 0)",
                        (os.path.join(hf.resource_path, "psi_pred"),))

    work_con.commit()

    worker.cpus = 2
    monkeypatch.setattr(launch_worker, "JOB_SIZE_COFACTOR", 2)
    monkeypatch.setattr(launch_worker.Worker, "prepare_all_by_all", lambda *_: [136, []])
    monkeypatch.setattr(launch_worker.Worker, "spawn_subjobs", lambda *_: [2, [], 1, 3])
    monkeypatch.setattr(launch_worker.Worker, "load_subjob", lambda *_: [4, []])

    def kill_worker(*args, **kwargs):
        worker.terminate("unit test kill")

    monkeypatch.setattr(br, "run_multicore_function", kill_worker)

    seqbuddy = hf.get_data("cteno_panxs")
    seqbuddy.write(os.path.join(worker.output, "foo.seqs"))

    alignment = hf.get_data("cteno_panxs_aln")
    alignment.write(os.path.join(worker.output, "foo.aln"))

    monkeypatch.setattr(Alb, "generate_msa", lambda *_, **__: alignment)

    # This first run is a full (large) job that will spawn subjobs
    with pytest.raises(SystemExit):
        worker.start()

    out, err = capsys.readouterr()
    assert "Terminating Worker_2 because of unit test kill" in out, print(out)

    # A second run is pulling a subjob off the queue
    work_cursor.execute("INSERT INTO "
                        "queue (hash, psi_pred_dir, master_id, align_m, align_p, trimal, gap_open, gap_extend) "
                        "VALUES ('2_3_foo', ?, 2, 'clustalo', '', 'gappyout 50 90 clean', 0, 0)",
                        (os.path.join(hf.resource_path, "psi_pred"),))
    work_con.commit()

    with pytest.raises(SystemExit):
        worker.start()

    out, err = capsys.readouterr()
    assert "Terminating Worker_3 because of unit test kill" in out, print(out)


def test_worker_check_master(hf, capsys):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    temp_dir.copy_to("%sheartbeat_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path)
    worker.last_heartbeat_from_master = time.time() + 100
    assert worker.check_masters(20) is None

    worker.last_heartbeat_from_master = time.time() - 200
    with pytest.raises(SystemExit):
        worker.check_masters(20)
    out, err = capsys.readouterr()
    assert "Terminating Worker_None because of 2 min, 0 sec of master inactivity (spent 20% time idle)" in out


def test_worker_masterclear(hf, capsys):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    temp_dir.copy_to("%sheartbeat_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path)
    masterclear = temp_dir.subfile("MasterClear")
    assert os.path.isfile(worker.masterclear_path), print(worker.masterclear_path)

    with open(masterclear, "w") as ofile:
        ofile.write("12345")
    assert worker.masterclear() == 12345
    launch_worker.sys.stdout.flush()
    out, err = capsys.readouterr()
    assert out == "\r\rMasterClear signal 12345\n"
    assert not os.path.isfile(masterclear)

    with open(masterclear, "w") as ofile:
        ofile.write("Foo")
    assert worker.masterclear() == worker.max_wait
    out, err = capsys.readouterr()
    assert out == "\r\rinvalid literal for int() with base 10: 'Foo'\n"


def test_worker_clean_dead_threads(hf):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    temp_dir.copy_to("%sheartbeat_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path)

    hb_con = sqlite3.connect(os.path.join(temp_dir.path, "heartbeat_db.sqlite"))
    hb_cursor = hb_con.cursor()
    hb_cursor.execute("INSERT INTO heartbeat (thread_type, pulse) VALUES ('master', %s)" % round(time.time()))
    master_id = hb_cursor.lastrowid
    hb_cursor.execute("INSERT INTO heartbeat (thread_type, pulse) VALUES ('worker', %s)" % round(time.time()))
    worker_id = hb_cursor.lastrowid
    hb_con.commit()

    work_con = sqlite3.connect(os.path.join(temp_dir.path, "work_db.sqlite"))
    work_cursor = work_con.cursor()
    work_cursor.execute("INSERT INTO "
                        "waiting (hash, master_id) "
                        "VALUES ('foo', %s)" % master_id)
    work_cursor.execute("INSERT INTO "
                        "processing (hash, worker_id, master_id) "
                        "VALUES ('foo', ?, ?)", (worker_id, master_id,))
    # Orphan
    work_cursor.execute("INSERT INTO "
                        "complete (hash, worker_id, master_id) "
                        "VALUES ('baz', 101, 100)")
    work_con.commit()

    worker.clean_dead_threads(100)
    assert hb_cursor.execute("SELECT * FROM heartbeat WHERE thread_id=%s" % master_id).fetchone()
    assert work_cursor.execute("SELECT * FROM waiting WHERE master_id=%s" % master_id).fetchone()
    assert work_cursor.execute("SELECT * FROM processing WHERE worker_id=%s" % worker_id).fetchone()

    # No one waiting on complete
    work_cursor.execute("INSERT INTO "
                        "complete (hash, worker_id, master_id) "
                        "VALUES ('bar', ?, ?)", (worker_id, master_id,))
    work_con.commit()

    worker.clean_dead_threads(-100)
    assert not hb_cursor.execute("SELECT * FROM heartbeat WHERE thread_id=%s" % master_id).fetchone()
    assert not work_cursor.execute("SELECT * FROM waiting WHERE master_id=%s" % master_id).fetchone()
    assert not work_cursor.execute("SELECT * FROM processing WHERE worker_id=%s" % worker_id).fetchone()


def test_worker_fetch_queue_job(hf):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path)

    work_con = sqlite3.connect(os.path.join(temp_dir.path, "work_db.sqlite"))
    work_cursor = work_con.cursor()
    work_cursor.execute("INSERT INTO "
                        "queue (hash, psi_pred_dir, master_id, align_m, align_p, trimal, gap_open, gap_extend) "
                        "VALUES ('foo', './', 2, '', '', 'gappyout 50 90 clean', 0, 0)")
    work_con.commit()
    queued_job = worker.fetch_queue_job()
    assert queued_job == ['foo', './', 2, '', '', ['gappyout', 50.0, 90.0, 'clean'], 0, 0]
    assert not work_cursor.execute("SELECT * FROM queue").fetchall()


def test_worker_prepare_psipred_dfs(hf, capsys):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path)
    seqbuddy = hf.get_data("cteno_panxs")
    psipred_dfs = worker.prepare_psipred_dfs(seqbuddy, os.path.join(hf.resource_path, "psi_pred"))
    assert len(seqbuddy) == len(psipred_dfs)
    assert str(psipred_dfs["Hvu-PanxβO"].iloc[0]) == """\
indx              1
aa                M
ss                C
coil_prob     0.999
helix_prob    0.001
sheet_prob    0.001
Name: 0, dtype: object"""

    capsys.readouterr()
    seqbuddy.records[0].id = "Foo"
    with pytest.raises(SystemExit):
        worker.prepare_psipred_dfs(seqbuddy, os.path.join(hf.resource_path, "psi_pred"))
    out, err = capsys.readouterr()
    assert "Terminating Worker_None because of missing psi ss2 file" in out


def test_worker_update_psipred(hf):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path)

    align = hf.get_data("cteno_panxs_aln")
    align = Alb.extract_regions(align, "105:115")
    align = Alb.pull_records(align, "Bfo-PanxαA|Bfr-PanxαD")
    ss2_dfs = hf.get_data("ss2_dfs")
    ss2_dfs = {"Bfo-PanxαA": ss2_dfs["Bfo-PanxαA"], "Bfr-PanxαD": ss2_dfs["Bfr-PanxαD"]}
    ss2_dfs["Bfo-PanxαA"] = ss2_dfs["Bfo-PanxαA"].iloc[47:56]
    for indx, new in [(47, 1), (48, 2), (49, 3), (50, 4), (51, 5), (52, 6), (53, 7), (54, 8), (55, 9)]:
        ss2_dfs["Bfo-PanxαA"].set_value(indx, "indx", new)
    ss2_dfs["Bfo-PanxαA"] = ss2_dfs["Bfo-PanxαA"].reset_index(drop=True)

    ss2_dfs["Bfr-PanxαD"] = ss2_dfs["Bfr-PanxαD"].iloc[44:53]
    for indx, new in [(44, 1), (45, 2), (46, 3), (47, 4), (48, 5), (49, 6), (50, 7), (51, 8), (52, 9)]:
        ss2_dfs["Bfr-PanxαD"].set_value(indx, "indx", new)

    ss2_dfs["Bfr-PanxαD"] = ss2_dfs["Bfr-PanxαD"].reset_index(drop=True)

    ss2_dfs = worker.update_psipred(align, ss2_dfs, "msa")
    assert str(ss2_dfs["Bfo-PanxαA"]) == """\
   indx aa ss  coil_prob  helix_prob  sheet_prob
0     0  S  H      0.034       0.966       0.003
1     1  Q  H      0.071       0.926       0.004
2     2  M  H      0.371       0.649       0.003
3     3  W  C      0.802       0.211       0.004
4     4  S  C      0.852       0.151       0.010
5     5  Q  C      0.765       0.253       0.009
6     8  D  C      0.733       0.283       0.011
7     9  D  C      0.890       0.126       0.014
8    10  A  C      0.914       0.085       0.028"""

    align = Alb.trimal(align, "all")
    ss2_dfs = worker.update_psipred(align, ss2_dfs, "trimal")
    assert str(ss2_dfs["Bfo-PanxαA"]) == """\
   indx aa ss  coil_prob  helix_prob  sheet_prob
0     0  S  H      0.034       0.966       0.003
1     3  W  C      0.802       0.211       0.004
2     4  S  C      0.852       0.151       0.010
3     5  Q  C      0.765       0.253       0.009
4     8  D  C      0.733       0.283       0.011
5     9  D  C      0.890       0.126       0.014
6    10  A  C      0.914       0.085       0.028""", print(str(ss2_dfs["Bfo-PanxαA"]))

    with pytest.raises(ValueError) as err:
        worker.update_psipred(align, ss2_dfs, "foo")

    assert "Unrecognized mode 'foo': select from ['msa', 'trimal']" in str(err)


def test_worker_trimal():
    temp_dir = br.TempDir()
    worker = launch_worker.Worker(temp_dir.path)
    align = Alb.AlignBuddy("""\
>A
MSTGTC-------
>B
M---TC-------
>C
M---TC---AILP
>D
-STP---YWAILP
""", in_format="fasta")

    seqbuddy = Sb.SeqBuddy(Alb.make_copy(align).records(), in_format="fasta")
    seqbuddy = Sb.clean_seq(seqbuddy)

    # Don't modify if any sequence is reduced to nothing
    trimal = worker.trimal(seqbuddy, [0.3], Alb.make_copy(align))
    assert str(trimal) == str(align)

    align = Alb.AlignBuddy("""\
>A
MSTGTC-------
>B
M---TC-------
>C
M---TC---AILP
>D
-STPTC-YWAILP
""", in_format="fasta")

    # Don't modify if average sequence length is reduced by more than half
    trimal = worker.trimal(seqbuddy, [0.3], Alb.make_copy(align))
    assert str(trimal) == str(align)

    # Remove some gaps
    trimal = worker.trimal(seqbuddy, ["all", 0.3, 0.55, "clean"], Alb.make_copy(align))
    assert str(trimal) == """\
>A
MSTGTC----
>B
M---TC----
>C
M---TCAILP
>D
-STPTCAILP
"""


def test_worker_prepare_all_by_all(hf, monkeypatch):
    temp_dir = br.TempDir()
    monkeypatch.setattr(br, "cpu_count", lambda *_: 24)
    worker = launch_worker.Worker(temp_dir.path)
    seqbuddy = hf.get_data("cteno_panxs")
    seqbuddy = Sb.pull_recs(seqbuddy, "Oma")  # Only 4 records, which means 6 comparisons
    ss2_dfs = hf.get_data("ss2_dfs")

    data_len, data = worker.prepare_all_by_all(seqbuddy, ss2_dfs)
    assert data_len == 6
    assert len(data) == 6
    assert data[0][0][0:2] == ('Oma-PanxαA', 'Oma-PanxαB')

    seqbuddy = hf.get_data("cteno_panxs")  # 134 records = 8911 comparisons
    data_len, data = worker.prepare_all_by_all(seqbuddy, ss2_dfs)
    assert data_len == 8911
    assert len(data[0]) == int(ceil(data_len / (br.cpu_count() - 1)))


def test_worker_process_final_results(hf, monkeypatch):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)

    work_con = sqlite3.connect(os.path.join(temp_dir.path, "work_db.sqlite"))
    work_cursor = work_con.cursor()

    worker = launch_worker.Worker(temp_dir.path)
    aln_file = temp_dir.subfile(os.path.join(worker.output, "foo.aln"))
    seqs_file = temp_dir.subfile(os.path.join(worker.output, "foo.seqs"))
    graph_file = os.path.join(worker.output, "foo.graph")
    subjob_dir = temp_dir.subdir(os.path.join(worker.output, "foo"))

    # Use worker id = 2 and master id = 3
    worker.data_file = os.path.join(temp_dir.path, ".Worker_2.dat")

    work_cursor.execute("INSERT INTO "
                        "waiting (hash, master_id) "
                        "VALUES ('foo', 3)")
    work_cursor.execute("INSERT INTO "
                        "processing (hash, worker_id, master_id) "
                        "VALUES ('foo', 2, 3)")
    work_cursor.execute("INSERT INTO "
                        "processing (hash, worker_id, master_id) "
                        "VALUES ('1_2_foo', 2, 3)")
    work_cursor.execute("INSERT INTO "
                        "complete (hash, worker_id, master_id) "
                        "VALUES ('1_2_foo', 2, 3)")
    work_con.commit()

    # Start by processing an empty result
    with open(worker.data_file, "w") as ifile:
        ifile.write("seq1,seq2,subsmat,psi")

    assert worker.process_final_results("foo", 3, 1, 1) is None
    assert work_cursor.execute("SELECT * FROM waiting WHERE hash='foo'").fetchone()
    assert work_cursor.execute("SELECT * FROM processing WHERE hash='foo'").fetchone()
    assert work_cursor.execute("SELECT * FROM processing WHERE hash LIKE '%%_foo'").fetchone()
    assert not work_cursor.execute("SELECT * FROM complete WHERE hash='foo'").fetchone()
    assert work_cursor.execute("SELECT * FROM complete WHERE hash LIKE '%%_foo'").fetchone()
    assert os.path.isfile(aln_file)
    assert os.path.isfile(seqs_file)
    assert not os.path.isfile(graph_file)
    assert os.path.isdir(subjob_dir)

    # Confirm that the sim_scores will be collected from process_subjob if num_subjobs > 1
    sim_scores = launch_worker.pd.read_csv(worker.data_file, index_col=False)
    monkeypatch.setattr(launch_worker.Worker, "process_subjob", lambda *_: sim_scores)

    assert worker.process_final_results("foo", 3, 1, 4) is None
    assert work_cursor.execute("SELECT * FROM waiting WHERE hash='foo'").fetchone()
    assert work_cursor.execute("SELECT * FROM processing WHERE hash='foo'").fetchone()
    assert not work_cursor.execute("SELECT * FROM complete WHERE hash='foo'").fetchone()
    assert os.path.isfile(aln_file)
    assert os.path.isfile(seqs_file)
    assert not os.path.isfile(graph_file)
    assert os.path.isdir(subjob_dir)

    # Process an actual result
    result = """
Oma-PanxαA,Oma-PanxαB,0.24409184734877187,0.48215666893146514
Oma-PanxαA,Oma-PanxαC,0.5135617078978646,0.7301561315022769
Oma-PanxαA,Oma-PanxαD,0.587799312776859,0.7428144752048478
Oma-PanxαB,Oma-PanxαC,0.23022898459442326,0.5027489193831555
Oma-PanxαB,Oma-PanxαD,0.2382962711779895,0.44814103743455824
Oma-PanxαC,Oma-PanxαD,0.47123287364498306,0.6647632735397735
"""
    with open(worker.data_file, "a") as ifile:
        ifile.write(result)
    assert worker.process_final_results("foo", 3, 1, 1) is None
    assert work_cursor.execute("SELECT * FROM waiting WHERE hash='foo'").fetchone()
    assert not work_cursor.execute("SELECT * FROM processing WHERE hash='foo'").fetchone()
    assert work_cursor.execute("SELECT * FROM complete WHERE hash='foo'").fetchone()
    assert not work_cursor.execute("SELECT * FROM complete WHERE hash LIKE '%%_foo'").fetchone()
    assert os.path.isfile(aln_file)
    assert os.path.isfile(seqs_file)
    assert os.path.isfile(graph_file)
    assert not os.path.isdir(subjob_dir)

    # Process a result with no masters waiting around
    work_cursor.execute("DELETE FROM waiting WHERE hash='foo'")
    work_cursor.execute("DELETE FROM complete WHERE hash='foo'")
    work_con.commit()
    os.remove(seqs_file)
    assert worker.process_final_results("foo", 3, 1, 1) is None
    assert not work_cursor.execute("SELECT * FROM waiting WHERE hash='foo'").fetchone()
    assert not work_cursor.execute("SELECT * FROM processing WHERE hash='foo'").fetchone()
    assert not work_cursor.execute("SELECT * FROM processing WHERE hash LIKE '%%_foo'").fetchone()
    assert not work_cursor.execute("SELECT * FROM complete WHERE hash='foo'").fetchone()
    assert not work_cursor.execute("SELECT * FROM complete WHERE hash LIKE '%%_foo'").fetchone()
    assert not os.path.isfile(aln_file)
    assert not os.path.isfile(seqs_file)
    assert not os.path.isfile(graph_file)


def test_worker_spawn_subjobs(hf, monkeypatch):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)

    work_con = sqlite3.connect(os.path.join(temp_dir.path, "work_db.sqlite"))
    work_cursor = work_con.cursor()

    worker = launch_worker.Worker(temp_dir.path)
    worker.heartbeat.id = 1
    subjob_dir = os.path.join(worker.output, "foo")

    # Set max job size at 4
    worker.cpus = 2
    monkeypatch.setattr(launch_worker, "JOB_SIZE_COFACTOR", 2)

    ss2_dfs = hf.get_data("ss2_dfs")

    pairs = [('Bch-PanxαA', 'Bch-PanxαB'), ('Bch-PanxαA', 'Bch-PanxαC'), ('Bch-PanxαA', 'Bch-PanxαD'),
             ('Bch-PanxαA', 'Bch-PanxαE'), ('Bch-PanxαB', 'Bch-PanxαC'), ('Bch-PanxαB', 'Bch-PanxαD'),
             ('Bch-PanxαB', 'Bch-PanxαE'), ('Bch-PanxαC', 'Bch-PanxαD'), ('Bch-PanxαC', 'Bch-PanxαE'),
             ('Bch-PanxαD', 'Bch-PanxαE')]

    for indx, pair in enumerate(pairs):
        pairs[indx] = (pair[0], pair[1], ss2_dfs[pair[0]], ss2_dfs[pair[1]])

    data = [pairs[i:i + 5] for i in range(0, len(pairs), 5)]  # This gives two groups of five

    data_len, data, subjob_num, num_subjobs = worker.spawn_subjobs("foo", data, ss2_dfs, 3, -5, 0)

    assert data_len == len(data[0]) == 2
    assert len(data[1]) == 2
    assert len(data) == 2
    assert subjob_num == 1
    assert num_subjobs == 3

    assert os.path.isfile(os.path.join(subjob_dir, "Bch-PanxαA.ss2"))
    with open(os.path.join(subjob_dir, "1_of_3.txt"), "r") as ifile:
        assert ifile.read() == """\
Bch-PanxαA Bch-PanxαB
Bch-PanxαA Bch-PanxαC
Bch-PanxαA Bch-PanxαD
Bch-PanxαA Bch-PanxαE
"""

    with open(os.path.join(subjob_dir, "3_of_3.txt"), "r") as ifile:
        assert ifile.read() == """\
Bch-PanxαC Bch-PanxαE
Bch-PanxαD Bch-PanxαE
"""

    queue = work_cursor.execute("SELECT * FROM queue").fetchall()
    assert len(queue) == 2
    hash_id, psi_pred_dir, master_id, align_m, align_p, trimal, gap_open, gap_extend = queue[0]
    assert hash_id == "2_3_foo"
    assert psi_pred_dir == subjob_dir
    assert master_id == 3
    assert align_m is align_p is trimal is None
    assert gap_open, gap_extend == (-5, 0)

    processing = work_cursor.execute("SELECT * FROM processing").fetchall()
    assert len(processing) == 1

    hash_id, worker_id, master_id = processing[0]
    assert hash_id == "1_3_foo"
    assert worker_id == worker.heartbeat.id
    assert master_id == 3


def test_worker_load_subjob(hf):
    temp_dir = br.TempDir()
    worker = launch_worker.Worker(temp_dir.path)
    temp_dir.subdir(".worker_output/foo")
    worker.cpus = 2
    ss2_dfs = hf.get_data("ss2_dfs")
    subjob_dir = os.path.join(worker.output, "foo")

    with open(os.path.join(subjob_dir, "2_of_3.txt"), "w") as ofile:
        ofile.write("""\
Bch-PanxαA Bch-PanxαB
Bch-PanxαA Bch-PanxαC
Bch-PanxαA Bch-PanxαD
Bch-PanxαA Bch-PanxαE
""")

    data_len, data = worker.load_subjob("foo", 2, 3, ss2_dfs)
    assert data_len == 4
    assert len(data) == 2


def test_worker_process_subjob(hf):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path)
    worker.heartbeat.id = 1

    subjob_dir = temp_dir.subdir(".worker_output/foo")

    work_con = sqlite3.connect(os.path.join(temp_dir.path, "work_db.sqlite"))
    work_cursor = work_con.cursor()
    work_cursor.execute("INSERT INTO "
                        "processing (hash, worker_id, master_id) "
                        "VALUES ('2_3_foo', 1, 3)")
    work_con.commit()

    # First test what happens when there are remaining subjobs
    sim_scores = """\
seq1,seq2,subsmat,psi
Oma-PanxαA,Oma-PanxαB,0.2440918473487719,0.4821566689314651
Oma-PanxαA,Oma-PanxαC,0.5135617078978646,0.7301561315022769
"""
    sim_scores_df = pd.read_csv(StringIO(sim_scores))
    sim_scores_df = worker.process_subjob("foo", sim_scores_df, 2, 3, 3)
    assert sim_scores_df.empty
    assert work_cursor.execute("SELECT * FROM complete").fetchone()[0] == "2_3_foo"
    assert not work_cursor.execute("SELECT * FROM processing").fetchall()

    with open(os.path.join(subjob_dir, "2_of_3.sim_df"), "r") as ifile:
        assert ifile.read() == sim_scores

    # Now test the final processing after all subjobs complete
    work_cursor.execute("INSERT INTO "
                        "complete (hash, worker_id, master_id) "
                        "VALUES ('1_3_foo', 1, 3)")

    work_cursor.execute("INSERT INTO "
                        "processing (hash, worker_id, master_id) "
                        "VALUES ('3_3_foo', 1, 3)")
    work_con.commit()

    with open(os.path.join(subjob_dir, "1_of_3.sim_df"), "w") as ofile:
        ofile.write("""\
seq1,seq2,subsmat,psi
Oma-PanxαA,Oma-PanxαD,0.587799312776859,0.7428144752048478
Oma-PanxαB,Oma-PanxαC,0.2302289845944233,0.5027489193831555
""")

    sim_scores = """\
seq1,seq2,subsmat,psi
Oma-PanxαB,Oma-PanxαD,0.2382962711779895,0.44814103743455824
Oma-PanxαC,Oma-PanxαD,0.4712328736449831,0.6647632735397735
"""
    sim_scores_df = pd.read_csv(StringIO(sim_scores))
    sim_scores_df = worker.process_subjob("foo", sim_scores_df, 3, 3, 3)

    assert str(sim_scores_df) == """\
         seq1        seq2   subsmat       psi
0  Oma-PanxαA  Oma-PanxαD  0.587799  0.742814
1  Oma-PanxαB  Oma-PanxαC  0.230229  0.502749
0  Oma-PanxαA  Oma-PanxαB  0.244092  0.482157
1  Oma-PanxαA  Oma-PanxαC  0.513562  0.730156
0  Oma-PanxαB  Oma-PanxαD  0.238296  0.448141
1  Oma-PanxαC  Oma-PanxαD  0.471233  0.664763"""


def test_worker_terminate(hf, monkeypatch, capsys):
    monkeypatch.setattr(launch_worker.rdmcl.HeartBeat, "end", lambda *_: True)
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    work_con = sqlite3.connect(os.path.join(temp_dir.path, "work_db.sqlite"))
    work_cursor = work_con.cursor()
    work_cursor.execute("INSERT INTO "
                        "processing (hash, worker_id, master_id) "
                        "VALUES ('2_3_foo', 1, 3)")
    work_cursor.execute("INSERT INTO "
                        "processing (hash, worker_id, master_id) "
                        "VALUES ('1_3_foo', 2, 3)")
    work_con.commit()

    worker = launch_worker.Worker(temp_dir.path)
    worker.heartbeat.id = 1
    worker.data_file = temp_dir.subfile(".Worker_1.dat")
    assert os.path.isfile(worker.data_file)
    worker.worker_file = temp_dir.subfile("Worker_1")
    assert os.path.isfile(worker.worker_file)

    with pytest.raises(SystemExit):
        worker.terminate("unit test signal")
    out, err = capsys.readouterr()
    assert "Terminating Worker_1 because of unit test signal" in out
    assert not os.path.isfile(worker.data_file)
    assert not os.path.isfile(worker.worker_file)
    assert not work_cursor.execute("SELECT * FROM processing WHERE worker_id=1").fetchall()
    assert work_cursor.execute("SELECT * FROM processing WHERE worker_id=2").fetchall()


def test_score_sequences(hf):
    outfile = br.TempFile()
    alb_obj = hf.get_data("cteno_panxs_aln")
    # Grab subset of alignment for manual calculation
    # Bfo-PanxαA   SQMWSQ--DDA
    # Bfr-PanxαD   V--RQIVVGGP
    alb_obj = Alb.extract_regions(alb_obj, "105:115")

    ss2_dfs = hf.get_data("ss2_dfs")
    ss2_dfs = {"Bfo-PanxαA": ss2_dfs["Bfo-PanxαA"], "Bfr-PanxαD": ss2_dfs["Bfr-PanxαD"]}
    # Update the ss2 dfs according to the alignment subsequences extracted
    '''
    indx aa ss  coil_prob  helix_prob  sheet_prob
47     1  S  H      0.034       0.966       0.003
48     2  Q  H      0.071       0.926       0.004
49     3  M  H      0.371       0.649       0.003
50     4  W  C      0.802       0.211       0.004
51     5  S  C      0.852       0.151       0.010
52     6  Q  C      0.765       0.253       0.009
53     9  D  C      0.733       0.283       0.011
54    10  D  C      0.890       0.126       0.014
55    11  A  C      0.914       0.085       0.028
    '''
    '''
    indx aa ss  coil_prob  helix_prob  sheet_prob
44     1  V  E      0.136       0.196       0.556
45     4  R  E      0.178       0.357       0.553
46     5  Q  H      0.157       0.530       0.525
47     6  I  H      0.114       0.771       0.319
48     7  V  H      0.217       0.675       0.212
49     8  V  C      0.461       0.443       0.166
50     9  G  C      0.837       0.040       0.077
51    10  G  C      0.940       0.008       0.063
52    11  P  C      0.606       0.015       0.402
    '''

    ss2_dfs["Bfo-PanxαA"] = ss2_dfs["Bfo-PanxαA"].iloc[47:56]
    for indx, new in [(47, 1), (48, 2), (49, 3), (50, 4), (51, 5), (52, 6), (53, 9), (54, 10), (55, 11)]:
        ss2_dfs["Bfo-PanxαA"].set_value(indx, "indx", new)

    ss2_dfs["Bfr-PanxαD"] = ss2_dfs["Bfr-PanxαD"].iloc[44:53]
    for indx, new in [(44, 1), (45, 4), (46, 5), (47, 6), (48, 7), (49, 8), (50, 9), (51, 10), (52, 11)]:
        ss2_dfs["Bfr-PanxαD"].set_value(indx, "indx", new)

    gap_open = -5
    gap_extend = 0

    # For score, subsmat = -0.363
    launch_worker.score_sequences([("Bfo-PanxαA", "Bfr-PanxαD", ss2_dfs["Bfo-PanxαA"], ss2_dfs["Bfr-PanxαD"])],
                                  [alb_obj, gap_open, gap_extend, outfile.path])

    assert outfile.read() == "\nBfo-PanxαA,Bfr-PanxαD,-0.3627272727272728,0.4183636363636363"


# #########  User Interface  ########## #
parser = argparse.ArgumentParser(prog="launch_worker", description="",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-wdb", "--workdb", action="store", default=os.getcwd(),
                    help="Specify the directory where sqlite databases will be fed by RD-MCL", )
parser.add_argument("-hr", "--heart_rate", type=int, default=60,
                    help="Specify how often the worker should check in")
parser.add_argument("-mw", "--max_wait", action="store", type=int, default=120,
                    help="Specify the maximum time a worker will stay alive without seeing a master")
parser.add_argument("-log", "--log", help="Stream log data one line at a time", action="store_true")
parser.add_argument("-q", "--quiet", help="Suppress all output", action="store_true")

in_args = parser.parse_args([])


def test_argparse_init(monkeypatch):
    out_dir = br.TempDir()
    argv = ['launch_worker.py', '--workdb', out_dir.path]
    monkeypatch.setattr(launch_worker.sys, "argv", argv)
    temp_in_args = launch_worker.argparse_init()
    assert temp_in_args.workdb == out_dir.path
    assert temp_in_args.heart_rate == 60
    assert temp_in_args.max_wait == 120
    assert not temp_in_args.log
    assert not temp_in_args.quiet


def test_main(monkeypatch, capsys):
    out_dir = br.TempDir()
    argv = ['launch_worker.py', '--workdb', out_dir.path, "--max_wait", "1"]
    monkeypatch.setattr(launch_worker.sys, "argv", argv)

    with pytest.raises(SystemExit):
        launch_worker.main()
    out, err = capsys.readouterr()
    assert 'Terminating Worker_1 because of 1 sec of master inactivity' in out

    workdb_con = sqlite3.connect(os.path.join(out_dir.path, "work_db.sqlite"))
    workdb_cursor = workdb_con.cursor()

    tables = {'queue': [(0, 'hash', 'TEXT', 0, None, 1),
                        (1, 'psi_pred_dir', 'TEXT', 0, None, 0),
                        (2, 'master_id', 'INTEGER', 0, None, 0),
                        (3, 'align_m', 'TEXT', 0, None, 0),
                        (4, 'align_p', 'TEXT', 0, None, 0),
                        (5, 'trimal', 'TEXT', 0, None, 0),
                        (6, 'gap_open', 'FLOAT', 0, None, 0),
                        (7, 'gap_extend', 'FLOAT', 0, None, 0)],
              'processing': [(0, 'hash', 'TEXT', 0, None, 1),
                             (1, 'worker_id', 'INTEGER', 0, None, 0),
                             (2, 'master_id', 'INTEGER', 0, None, 0)],
              'complete': [(0, 'hash', 'TEXT', 0, None, 1),
                           (1, 'worker_id', 'INTEGER', 0, None, 0),
                           (2, 'master_id', 'INTEGER', 0, None, 0)],
              'waiting': [(0, 'hash', 'TEXT', 0, None, 0),
                          (1, 'master_id', 'INTEGER', 0, None, 0)]}

    workdb_tables = workdb_cursor.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
    for table in tables:
        assert (table,) in workdb_tables

    for table, fields in tables.items():
        fields_query = workdb_cursor.execute("PRAGMA table_info(%s)" % table).fetchall()
        assert fields_query == fields

    hb_db_con = sqlite3.connect(os.path.join(out_dir.path, "heartbeat_db.sqlite"))
    hb_db_cursor = hb_db_con.cursor()
    tables = {'heartbeat': [(0, 'thread_id', 'INTEGER', 0, None, 1),
                            (1, 'thread_type', 'TEXT', 0, None, 0),
                            (2, 'pulse', 'INTEGER', 0, None, 0)]}

    hb_db_tables = hb_db_cursor.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
    for table in tables:
        assert (table,) in hb_db_tables

    for table, fields in tables.items():
        fields_query = hb_db_cursor.execute("PRAGMA table_info(%s)" % table).fetchall()
        assert fields_query == fields

    # Test quiet
    argv = ['launch_worker.py', '--workdb', out_dir.path, "--max_wait", "1", "--quiet"]
    monkeypatch.setattr(launch_worker.sys, "argv", argv)
    capsys.readouterr()
    with pytest.raises(SystemExit):
        launch_worker.main()
    out, err = capsys.readouterr()
    assert not (out + err)

    out, err = capsys.readouterr()
    assert not (out + err)

    # Test logging mode (line breaks are inserted between start and termination messages)
    argv = ['launch_worker.py', '--workdb', out_dir.path, "--max_wait", "2", "--log"]
    monkeypatch.setattr(launch_worker.sys, "argv", argv)
    capsys.readouterr()
    with pytest.raises(SystemExit):
        launch_worker.main()

    out, err = capsys.readouterr()
    assert "Starting Worker_3\n" in out and "Terminating Worker_3 because of 2 sec of master inactivity" in out

    # Test termination types
    monkeypatch.setattr(launch_worker.helpers, "dummy_func", mock_valueerror)
    argv = ['launch_worker.py', '--workdb', out_dir.path, "--max_wait", "1"]
    monkeypatch.setattr(launch_worker.sys, "argv", argv)

    with pytest.raises(SystemExit):
        launch_worker.main()
    out, err = capsys.readouterr()
    assert 'Terminating Worker_7 because of too many Worker crashes' in out

    monkeypatch.setattr(launch_worker.helpers, "dummy_func", mock_keyboardinterupt)

    with pytest.raises(SystemExit):
        launch_worker.main()
    out, err = capsys.readouterr()
    assert 'Terminating Worker_8 because of KeyboardInterrupt' in out
