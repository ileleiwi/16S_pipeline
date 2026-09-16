# Memory Bank

This directory is the **persistent context** for the 16S Amplicon Sequencing
Pipeline project. It is designed to be read by Cline (and humans) at the start
of every task to rebuild full project context, since Cline's memory resets
between sessions.

## How to Use
1. **At the start of every task:** read ALL files in this directory.
2. **After any significant change:** update the relevant file(s).
3. **Keep it current:** stale context is worse than no context.

## Files

| File | Purpose |
|---|---|
| `projectbrief.md` | Foundational goals, scope, and success criteria |
| `productContext.md` | Why the project exists; problems it solves; domain background |
| `systemPatterns.md` | Architecture, key technical decisions, design patterns |
| `techContext.md` | Technologies, tools, setup, constraints |
| `activeContext.md` | Current work focus, recent changes, next steps |
| `progress.md` | What works, what's left, known issues |

## Hierarchy
```
projectbrief.md   (foundation - rarely changes)
      |
      +-- productContext.md   (why / domain)
      +-- systemPatterns.md   (how / architecture)
      +-- techContext.md      (with what / tools)
              |
              +-- activeContext.md  (now / focus)
              +-- progress.md       (status / roadmap)
```

## Maintenance Rules
- `projectbrief.md` changes only when scope or goals change.
- `activeContext.md` and `progress.md` are updated most frequently.
- When a decision is made (e.g., denoiser choice), record it in
  `systemPatterns.md` and remove it from the open questions in
  `activeContext.md`.
- Never store secrets, credentials, or large data here.