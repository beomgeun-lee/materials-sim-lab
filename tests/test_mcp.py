"""MCP 서버 — 실제 MCP 클라이언트로 붙어 도구를 부른다 (D23)."""

from __future__ import annotations

import asyncio
import json
from typing import Any

import pytest
from mcp import Client

from msl import mcp_server

RECIPE = """
id: rcp-mcp-salt
name: 소금물
components:
  - {ref: "cas:7647-14-5", amount: 0.1 mol, state: aqueous}
  - {ref: "cas:7732-18-5", amount: 1 L, state: liquid}
conditions: {T: 25 °C}
"""


def call(tool: str, args: dict[str, Any] | None = None):
    async def go():
        async with Client(mcp_server.mcp) as c:
            return await c.call_tool(tool, args or {})

    res = asyncio.run(go())
    if res.is_error:
        raise ToolFailed(text(res))
    return res


class ToolFailed(Exception):
    pass


def text(res) -> str:
    return "".join(getattr(x, "text", "") for x in res.content)


def data(res):
    return res.structured_content.get("result", res.structured_content) if res.structured_content else json.loads(text(res))


@pytest.fixture(autouse=True)
def _tmp_dirs(tmp_path, monkeypatch):
    monkeypatch.setattr(mcp_server, "REPORTS", tmp_path / "reports")
    monkeypatch.setattr(mcp_server.wb, "USER_RECIPES", tmp_path / "recipes")


def test_lists_minimal_tools() -> None:
    async def go():
        async with Client(mcp_server.mcp) as c:
            return {t.name: t for t in (await c.list_tools()).tools}

    tools = asyncio.run(go())
    assert {"search_substance", "check_recipe", "run_recipe", "get_report", "list_reports", "list_assays",
            "predict_properties", "recommend", "l2_result", "suggest_next", "add_measurement", "list_campaigns"} <= set(tools)
    assert tools["get_report"].annotations.read_only_hint is True
    assert tools["save_recipe"].annotations.read_only_hint is False


def test_list_assays_has_all_eleven() -> None:
    codes = {a["code"] for a in data(call("list_assays"))}
    assert codes == {"S0", *(f"A{i}" for i in range(1, 11))}


def test_search_substance_korean_alias() -> None:
    r = data(call("search_substance", {"text": "염산"}))
    assert r["formula"] == "HCl" and r["error"] is None


def test_check_recipe_reports_errors_in_korean() -> None:
    r = data(call("check_recipe", {"recipe_yaml": RECIPE.replace("0.1 mol", "-1 mol")}))
    assert r["ok"] is False and r["errors"]
    ok = data(call("check_recipe", {"recipe_yaml": RECIPE}))
    assert ok["ok"] and "A4" in {a["code"] for a in ok["assays"]}


def test_run_then_get_report() -> None:
    out = text(call("run_recipe", {"recipe_yaml": RECIPE}))
    first, _, body = out.partition("\n")
    rid = first.removeprefix("report_id: ")
    assert "pH" in body and "A4" in body
    assert [r["report_id"] for r in data(call("list_reports"))] == [rid]
    assert json.loads(text(call("get_report", {"report_id": rid, "fmt": "json"})))["recipe"]["id"] == "rcp-mcp-salt"


def test_run_example_file() -> None:
    assert "report_id:" in text(call("run_recipe", {"file": "calcite-water.yaml"}))


def test_report_id_rejects_paths() -> None:
    with pytest.raises(ToolFailed, match="리포트가 없음"):
        call("get_report", {"report_id": "../../.env"})


def test_save_recipe_refuses_example_id_and_needs_overwrite() -> None:
    ex = text(call("get_recipe", {"file": "calcite-water.yaml"}))
    with pytest.raises(ToolFailed, match="예제"):
        call("save_recipe", {"recipe_yaml": ex})
    assert data(call("save_recipe", {"recipe_yaml": RECIPE}))["file"] == "rcp-mcp-salt.yaml"
    with pytest.raises(ToolFailed):
        call("save_recipe", {"recipe_yaml": RECIPE})
    call("save_recipe", {"recipe_yaml": RECIPE, "overwrite": True})


def test_recommend_goal_errors_are_readable() -> None:
    with pytest.raises(ToolFailed, match="원소"):
        call("recommend", {"goal_yaml": "id: goal-x\nname: x\nrequired: [Xx]"})


def test_predict_properties() -> None:
    try:
        from msl.ml.predict import info

        info()
    except Exception:
        pytest.skip("L1 모델 없음 (msl ml train)")
    r = data(call("predict_properties", {"formulas": ["LiCoO2", "염산"]}))
    assert r[0]["formula"] == "LiCoO2" and "error" not in r[0]


def test_campaign_tools(tmp_path, monkeypatch) -> None:
    pytest.importorskip("baybe")
    from msl import suggest as sg

    monkeypatch.setattr(sg, "STORE", tmp_path / "campaigns")
    assert "li-mn-o-cathode.yaml" in {c["file"] for c in data(call("list_campaigns"))}
    with pytest.raises(ToolFailed, match="source"):
        call("add_measurement", {"campaign_file": "li-mn-o-cathode.yaml", "formula": "LiMnO2", "value": 0.0, "source": "추정"})
    with pytest.raises(ToolFailed, match="캠페인 파일이 없음"):
        call("suggest_next", {"campaign_file": "../../.env"})
    assert data(call("add_measurement", {"campaign_file": "li-mn-o-cathode.yaml", "formula": "LiMnO2", "value": 0.0,
                                         "source": "L2"}))["formula"] == "LiMnO2"
