import { useRef, useState } from 'react'


type Msg = { role: 'user'|'assistant', text: string }


export default function Chat() {
const [msgs, setMsgs] = useState<Msg[]>([])
const [input, setInput] = useState('')
const [yearMin, setYearMin] = useState<number | ''>('')
const [meshAny, setMeshAny] = useState('') // comma-separated
const ctrl = useRef<AbortController | null>(null)
const [busy, setBusy] = useState(false)


async function send() {
if (!input.trim()) return
setMsgs(m => [...m, { role: 'user', text: input }])


const body = {
message: input,
year_min: yearMin || undefined,
mesh_any: meshAny ? meshAny.split(',').map(s => s.trim()).filter(Boolean) : undefined,
}


setInput('')
setBusy(true)


ctrl.current?.abort()
ctrl.current = new AbortController()


const res = await fetch('/api/chat', {
method: 'POST',
headers: { 'Content-Type': 'application/json' },
body: JSON.stringify(body),
signal: ctrl.current.signal,
})


const reader = res.body?.getReader()
const decoder = new TextDecoder()
let assistant = ''


while (reader) {
const { value, done } = await reader.read()
if (done) break
assistant += decoder.decode(value, { stream: true })
setMsgs(m => {
const copy = [...m]
const last = copy[copy.length - 1]
if (!last || last.role !== 'assistant') copy.push({ role: 'assistant', text: assistant })
else copy[copy.length - 1] = { role: 'assistant', text: assistant }
return copy
})
}


setBusy(false)
}


return (
<div>
<div style={{ marginBottom: 12, display: 'flex', gap: 8 }}>
<input
placeholder="MeSH (comma-separated)"
value={meshAny}
onChange={e => setMeshAny(e.target.value)}
style={{ flex: 1 }}
/>
<input
placeholder="Year ≥"
value={yearMin}
onChange={e => setYearMin(e.target.value ? Number(e.target.value) : '')}
style={{ width: 100 }}
/>
</div>


<div style={{
border: '1px solid #ddd', borderRadius: 8, padding: 12, minHeight: 240,
whiteSpace: 'pre-wrap'
}}>
{msgs.map((m, i) => (
<div key={i} style={{ marginBottom: 10 }}>
<strong>{m.role === 'user' ? 'You' : 'Assistant'}:</strong> {m.text}
</div>
))}
</div>


<div style={{ display: 'flex', gap: 8, marginTop: 12 }}>
}